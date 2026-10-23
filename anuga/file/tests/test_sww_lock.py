"""Two runs writing one SWW file must fail loudly, not interleave (#232).

The writer re-opens the file in append mode at every yieldstep, so nothing at
the netCDF level stops a second run with the same name from appending its
frames between the first run's. The file then opens fine and every frame
holds plausible data; only the non-monotonic time variable gives it away.
Two guards: a sidecar lock taken when the file is created, and a refusal to
append a frame earlier than the last one on file.
"""
import os
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
import warnings

import anuga
from anuga import Domain, Reflective_boundary
from anuga.file.sww import (SWW_file, SWWFileInUseError, SWWTimeOrderError,
                            sww_lock_path)
from anuga.file.netcdf import NetCDFFile


def make_domain(datadir, name):
    points, vertices, boundary = anuga.rectangular_cross(4, 4, 1.0, 1.0)
    domain = Domain(points, vertices, boundary)
    domain.set_flow_algorithm('DE0')
    domain.set_datadir(datadir)
    domain.set_name(name)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 0.5)
    Br = Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain


class Test_sww_lock(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        for f in os.listdir(self.tmp):
            os.remove(os.path.join(self.tmp, f))
        os.rmdir(self.tmp)

    def test_lock_is_taken_on_create_and_released(self):
        domain = make_domain(self.tmp, 'locked')
        domain.initialise_storage()
        lock = sww_lock_path(domain.writer.filename)
        assert os.path.exists(lock), lock
        with open(lock) as f:
            assert 'pid=%d' % os.getpid() in f.read()
        domain.writer.release()
        assert not os.path.exists(lock)

    def test_a_live_process_holding_the_lock_is_refused(self):
        domain = make_domain(self.tmp, 'busy')
        lock = sww_lock_path(os.path.join(self.tmp, 'busy.sww'))
        # The parent (pytest's launcher, or the shell) is alive and is not us
        other = os.getppid()
        assert other != os.getpid()
        with open(lock, 'w') as f:
            f.write('pid=%d\nhost=%s\n' % (other, __import__('socket').gethostname()))
        with self.assertRaises(SWWFileInUseError) as cm:
            domain.initialise_storage()
        assert 'busy.sww' in str(cm.exception)
        assert 'process %d' % other in str(cm.exception)
        os.remove(lock)

    def test_a_stale_lock_is_taken_over_silently(self):
        domain = make_domain(self.tmp, 'stale')
        lock = sww_lock_path(os.path.join(self.tmp, 'stale.sww'))
        # A process that has certainly exited
        proc = subprocess.Popen([sys.executable, '-c', 'pass'])
        proc.wait()
        dead = proc.pid
        with open(lock, 'w') as f:
            f.write('pid=%d\nhost=%s\n' % (dead, __import__('socket').gethostname()))
        # Debris from a dead run must not warn: a suite running under
        # simplefilter('error') would fail on an unrelated stale lock.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            domain.initialise_storage()
        assert not [x for x in w if 'lock' in str(x.message)], \
            [str(x.message) for x in w]
        with open(lock) as f:
            assert 'pid=%d' % os.getpid() in f.read()

    def test_our_own_lock_is_taken_over_silently(self):
        """Re-running a model under the same name in one session is normal."""
        domain = make_domain(self.tmp, 'again')
        domain.initialise_storage()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            domain2 = make_domain(self.tmp, 'again')
            domain2.initialise_storage()
        assert not [x for x in w if 'lock' in str(x.message)], \
            [str(x.message) for x in w]
        assert os.path.exists(sww_lock_path(domain2.writer.filename))

    def test_a_frame_earlier_than_the_last_one_is_refused(self):
        domain = make_domain(self.tmp, 'backwards')
        for _ in domain.evolve(yieldstep=0.2, finaltime=0.4):
            pass
        domain.set_time(0.1)               # no frame at t = 0.1 to rewrite
        with self.assertRaises(SWWTimeOrderError) as cm:
            domain.store_timestep()
        assert 't=0.1' in str(cm.exception) and 't=0.4' in str(cm.exception)

    def test_rewriting_an_existing_frame_is_allowed(self):
        """A checkpoint resume, or a second evolve() re-storing its start,
        writes a time already on file: that frame is rewritten in place."""
        domain = make_domain(self.tmp, 'rewrite')
        for _ in domain.evolve(yieldstep=0.2, finaltime=0.4):
            pass
        domain.store_timestep()            # again at t = 0.4 (the last)
        domain.set_time(0.2)
        domain.store_timestep()            # again at t = 0.2 (an earlier one)
        with NetCDFFile(domain.writer.filename, 'r') as fid:
            t = [float(x) for x in fid.variables['time'][:]]
        assert t == [0.0, 0.2, 0.4], t

    def test_two_processes_cannot_write_one_file(self):
        """The real case: a second interpreter with the same name fails to
        start its output while the first is still running."""
        script = textwrap.dedent('''
            import os, sys, time
            import anuga
            from anuga import Domain, Reflective_boundary
            datadir, name, mode = sys.argv[1:4]
            points, vertices, boundary = anuga.rectangular_cross(4, 4, 1.0, 1.0)
            d = Domain(points, vertices, boundary)
            d.set_datadir(datadir); d.set_name(name)
            d.set_quantity('elevation', 0.0); d.set_quantity('stage', 0.5)
            Br = Reflective_boundary(d)
            d.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
            d.initialise_storage()
            if mode == 'hold':
                open(os.path.join(datadir, 'ready'), 'w').close()
                time.sleep(30)
        ''')
        holder = subprocess.Popen([sys.executable, '-c', script, self.tmp, 'shared', 'hold'],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            ready = os.path.join(self.tmp, 'ready')
            for _ in range(600):
                if os.path.exists(ready):
                    break
                time.sleep(0.1)
            assert os.path.exists(ready), 'holder never started'
            second = subprocess.run([sys.executable, '-c', script, self.tmp, 'shared', 'once'],
                                    capture_output=True, text=True)
            assert second.returncode != 0
            assert 'SWWFileInUseError' in second.stderr, second.stderr[-2000:]
            assert 'process %d' % holder.pid in second.stderr
        finally:
            holder.kill()
            holder.wait()

    def test_the_conflict_names_when_and_what_started_the_other_run(self):
        """A bare pid is little help when the other run is a notebook in
        another tab, so the message carries the start time and script the
        lock records."""
        domain = make_domain(self.tmp, 'named')
        lock = sww_lock_path(os.path.join(self.tmp, 'named.sww'))
        alive = os.getppid()                  # a pid that is certainly running
        with open(lock, 'w') as f:
            f.write('pid=%d\nhost=%s\nstarted=2026-09-28T09:15:02\n'
                    'script=/somewhere/ipykernel_launcher.py\n'
                    % (alive, __import__('socket').gethostname()))
        with self.assertRaises(SWWFileInUseError) as cm:
            domain.initialise_storage()
        msg = str(cm.exception)
        assert 'started 2026-09-28T09:15:02' in msg, msg
        assert 'ipykernel_launcher.py' in msg, msg
        os.remove(lock)

    def test_a_lock_without_those_fields_still_reports_the_process(self):
        """Locks written by an older ANUGA carry only pid and host."""
        domain = make_domain(self.tmp, 'bare')
        lock = sww_lock_path(os.path.join(self.tmp, 'bare.sww'))
        alive = os.getppid()
        with open(lock, 'w') as f:
            f.write('pid=%d\nhost=%s\n' % (alive, __import__('socket').gethostname()))
        with self.assertRaises(SWWFileInUseError) as cm:
            domain.initialise_storage()
        msg = str(cm.exception)
        assert 'process %d' % alive in msg, msg
        assert 'started' not in msg, msg
        os.remove(lock)


if __name__ == '__main__':
    unittest.main()

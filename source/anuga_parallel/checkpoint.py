"""
Procedures to support checkpointing

There is already checkpointing available in domain.

Setup with  domain.set_checkpointing(checkpoint_step, checkpoint_dir)

checkpoint_step: the number of yieldsteps between saving a checkpoint file
checkpoint_dir: the name of the directory where teh checkpoint files are stored. 


But if we are restarting a calculation there is no domain yet available, so we must 
read in the last stored domain. Do that via

domain = load_last_checkpoint_file(domain_name, checkpoint_dir)

"""

from anuga import send, receive, myid, numprocs



def load_checkpoint_file(domain_name = 'domain', checkpoint_dir = '.', time = None):
    
    from os.path import join

    if numprocs > 1:
        domain_name = domain_name+'_P{}_{}'.format(numprocs,myid)
            
    if time is None:
        # will pull out the last available time
        times = _get_checkpoint_times(domain_name, checkpoint_dir)

        times = list(times)
        times.sort()
        #print times
    else:
        times = [float(time)]
    
    for time in reversed(times):
        
        pickle_name = join(checkpoint_dir,domain_name)+'_'+str(time)+'.pickle'
        #print pickle_name
        
        try:
            import cPickle
            domain = cPickle.load(open(pickle_name, 'rb'))
            success = True
        except:
            success = False
            
        overall = success
        for cpu in range(numprocs):
            if myid != cpu:
                send(success,cpu)
                overall = overall and receive(cpu)  
        
        #print myid, overall, success
        
        if overall: break
    
    if not overall: raise Exception, "Unable to open checkpoint file"
    
    return domain


def _get_checkpoint_times(domain_name, checkpoint_dir):

    import os
    times = set()
    
    for (path, directory, filenames) in os.walk(checkpoint_dir):
        #print filenames
        #print directory
        
        if len(filenames) == 0:
            return None
        else:          
            for filename in filenames:
                filebase = os.path.splitext(filename)[0].rpartition("_")
                time = filebase[-1]
                domain_name_base = filebase[0]
                if domain_name_base == domain_name :
                    #print domain_name_base, time
                    times.add(float(time))


    #times.sort()
    #times = set(times)
    
    #print times
    combined = times
    for cpu in range(numprocs):
        if myid != cpu:
            send(times,cpu)
            rec = receive(cpu) 
            #print rec
            combined = combined & rec 
    
    #combined = list(combined).sort()
    #print combined
    
    return combined
    

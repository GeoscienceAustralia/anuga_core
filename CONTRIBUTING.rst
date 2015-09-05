
=================
Contributing code
=================

The preferred way to contribute to anuga is to fork the 
`main repository <http://github.com/GeoscienceAustralia/anuga_core>`_ on GitHub:

1. Fork the `project repository <http://github.com/GeoscienceAustralia/anuga_core>`_:
   click on the 'Fork' button near the top of the page. This creates
   a copy of the code under your account on the GitHub server.

2. Clone this copy to your local disk::

          $ git clone git@github.com:YourLogin/anuga_core.git
          $ cd anuga_core

3. Create a branch to hold your changes::

          $ git checkout -b my-feature

   and start making changes. Never work in the ``master`` branch!

4. Work on this copy on your computer using Git to do the version
   control. When you're done editing, do::

          $ git add modified_files
          $ git commit

   to record your changes in Git, then push them to GitHub with::

          $ git push -u origin my-feature

Finally, go to the web page of the your fork of the anuga_core repo,
and click 'Pull request' to send your changes to the maintainers for
review. request. This will send an email to the committers.

(If any of the above seems like magic to you, then look up the 
`Git documentation <http://git-scm.com/documentation>`_ on the web.)

It is recommended to check that your contribution complies with the
following rules before submitting a pull request:

-  All public methods should have informative docstrings with sample
   usage presented as doctests when appropriate.

-  All other tests pass when everything is rebuilt from scratch. On
   Unix-like systems, check with (from the toplevel source folder)::

          $ make

-  When adding additional functionality, provide at least one
   example script in the ``examples/`` folder. Have a look at other
   examples for reference. Examples should demonstrate why the new
   functionality is useful in practice.

-  At least one paragraph of narrative documentation with links to
   references in the literature (with PDF links when possible) and
   an example.

You can also check for common programming errors with the following
tools:

-  Code with good unittest coverage (at least 80%), check with::

          $ pip install nose coverage
          $ nosetests --with-coverage path/to/tests_for_package

-  No pyflakes warnings, check with::

           $ pip install pyflakes
           $ pyflakes path/to/module.py

-  No PEP8 warnings, check with::

           $ pip install pep8
           $ pep8 path/to/module.py

-  AutoPEP8 can help you fix some of the easy redundant errors::

           $ pip install autopep8
           $ autopep8 path/to/pep8.py

Bonus points for contributions that include a performance analysis with
a benchmark script and profiling output (please report on the mailing
list or on the GitHub issue).

Easy Issues
-----------

A great way to start contributing to anuga is to pick an item
from the list of `Easy issues <https://github.com/GeoscienceAustralia/anuga_core/issues?labels=Easy>`_
in the issue tracker. (Well there are none there yet, but we will be 
putting some up soon!) Resolving these issues allow you to start
contributing to the project without much prior knowledge. Your
assistance in this area will be greatly appreciated by the more
experienced developers as it helps free up their time to concentrate on
other issues.

Documentation
-------------

We are in the process of creating sphinx based documentation for anuga. 
Any help in setting this up will be gratefully accepted!

At present you will find the user_manual in the doc folder. 

We are glad to accept any sort of documentation: function docstrings,
reStructuredText documents (like this one), tutorials, etc.
reStructuredText documents live in the source code repository under the
doc/ directory.

When you are writing documentation, it is important to keep a good
compromise between mathematical and algorithmic details, and give
intuition to the reader on what the algorithm does. It is best to always
start with a small paragraph with a hand-waving explanation of what the
method does to the data and a figure (coming from an example)
illustrating it.


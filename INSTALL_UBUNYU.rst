= Installing current version of ANUGA research code on Ubuntu  14.04 =

The current development version of {{{anuga}}} now resides on {{{github}}} at 
[https://github.com/GeoscienceAustralia/anuga_core anuga_core]



== Install {{{git}}} ==

We use {{{git}}} to install {{{anuga}}} (subversion can still be used)

Use {{{apt-get}}} to install {{{git}}} via

{{{
sudo apt-get install git
}}}

== Install {{{anuga}}} ==

Now from your home directory, use subversion to download the {{{anuga}}} source via

{{{
git clone https://github.com/GeoscienceAustralia/anuga_core.git
}}}

We need to install a fairly large number of packages, but we have a script 
({{{tools/install_ubuntu.sh}}}) which will run through the installation process. 


=== Parallel Support ===
At this stage you can decide whether you want Parallel support or not. We support two versions of {{{MPI}}}, {{{mpich2}}} and {{{openmpi}}}

If you setup the environment variable  {{{ANUGA_PARALLEL}}} via:

{{{
export ANUGA_PARALLEL="mpich2"
}}}

or 
{{{
export ANUGA_PARALLEL="openmpi"
}}}

then the install script will load the  mpich2 or openmpi libraries and binaries respectively.

=== Running the installation script ===


Change into the newly downloaded {{{anuga_core}}} directory and run the installation script 
(this will take 5 to 10 minutes depending on your network connection)

{{{
cd anuga_core
bash tools/install_ubuntu.sh
}}}

If all has been successful then {{{anuga}}} should be installed.

== Check Installation ==
To check the installation run the over 1200 unit tests (which should take 5 - 10 minutes) via:

{{{
python runtests.py
}}}



Hopefully all the unit tests pass. As this is bleeding edge there are sometimes a small 
number of failures as this is a work in progress. 

== What Next ==

Have a look at the examples in the directory anuga_core/examples and the user_manual 
in the anuga_core/doc directory to see how to use anuga.

== Updating ==

From time to time you should update your version of anuga. This is fairly easy. 
From your {{{anuga_core}}} directory update the anuga code via the git command

{{{
git pull 
}}}

Then again from the {{{anuga_core}}} directory build and install the code and check the unit tests via

{{{
python setup.py build
sudo python setup.py install
python runtests.py
}}}

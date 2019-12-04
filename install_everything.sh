

# Install ANUGA as per docs

# Install SWIMM as per Stephen's instructions

git clone https://github.com/stoiver/Stormwater-Management-Model
cd Stormwater-Management-Model

# Switch Branch
git checkout toolkitapi

# Make shared library
cd build/Linux/
make libswmm5.so
make swmm5
sudo make install


# Install pyswmm

git clone  https://github.com/stoiver/pyswmm
cd pyswmm
git checkout linux
python setup.py install



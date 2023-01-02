# cpp-interface-osqp

Interface between DQ_Robotics and [OSQP](https://github.com/osqp/osqp)

## 1. Recommended installation process for OSQP

Check the official OSQP [instructions](https://osqp.org/docs/get_started/sources.html) 

Summary (Ubuntu):

```shell
cd ~/Downloads
git clone --recursive https://github.com/osqp/osqp
cd osqp
mkdir build && cd build
cmake -G "Unix Makefiles" ..
cmake --build .
cmake --build . --target install
```


## 2. Build the cpp-interface-osqp from sources (Tested only in MacOS ARM64)

```shell
cd ~/Downloads
git clone https://github.com/juanjqo/cpp-interface-osqp.git
cd cpp-interface-osqp
mkdir build
cd build
cmake ..
make -j16
sudo make install
```

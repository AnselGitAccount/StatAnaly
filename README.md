# StatAnaly: A Statisticaly Analysis Library

StatAnaly is a statisticaly analysis package written in C++. The header-only library is designed for rapid prototyping in numerical computations. 

StatAnaly is a work-in-progress in terms of features. Currently it supports the following operations:

1. Compute the sum of random variables of probability distributions.


## Quick Start Examples

There are three global static variables representing different sums:

- *cnvl* represents the simple sum, ie $R=X+Y$
- *cnvlSq* represents the sum of the squares, ie $R=X^2+Y^2$
- *cnvlSSqrt* represents the sum of the squares and then take square root, ie $R=\sqrt{X^2+Y^2}$

### Simple sum of probability distributions

The sum of two RVs of standard uniform distributions is a RV of Irwin-Hall distribution.

```cpp
disStdUniform su1{}, su2{};
probDistr* rsu = cnvl.go(su1, su2);
disIrwinHall* rih = dynamic_cast<disIrwinHall*>(rsu);
```

The sum of an array of RVs of standard uniform distributions is a RV of Irwin-Hall distribution.

```cpp
disStdUniform su1{}, su2{}, su3{}, su4{};
probDistr* rsu = convolve<disStdUniform>({su1, su2, su3, su4});
disIrwinHall* rih = dynamic_cast<disIrwinHall*>(rsu);
```

### Sum of the squares of probability distributions

The sum of two squares of RVs of normal distributions is a RV of Chi Square distribution.

```cpp
disNormal a(0,1);
disNormal b(0,1);
disChiSq s = *static_cast<disChiSq*>(cnvlSq.go(a,b));
```

### Sum of the Squares, and Then Take Square Root

The sum of two squares of RVs of normal distribution (with same variance) is a RV of Rician distribution.

```cpp
disNormal a(2,9);
disNormal b(3,9);
disRician s = *static_cast<disRician*>(cnvlSSqrt.go(a,b));
```

## Build & Install

### Requirements

The building and testing of the package require the following FOSS:

- CMake
- CTest
- Google Test

### Build StatAnaly

Make a build directory, and then compile the source codes.

```bash
    mkdir build
    cd build/
    cmake ../

    make
```

### Build a playground

A main function (in *playground.cpp*) is provided as a playground for experimentation. To build the executable:

```bash
    mkdir build
    cd build/
    cmake ../

    make playground
```

### Build & run unit tests

StatAnaly follows a test-driven development mindset. Unit tests are provided to ensure numerical accurarcy. StatAnaly supports CTest testing framework.

```bash
    mkdir build
    cd build/
    cmake ../

    make 
    make test
```

### Install StatAnaly

Install static library in */usr/local/lib/StatAnaly_x.x* and header files in */usr/local/include/StatAnaly_x.x*. Installation may require *sudo* privilege.

```bash
    mkdir build
    cd build/
    cmake ../

    sudo make install -j6
```

## Generate Doxygen

Source code is commented in a Doxygen recognizable style.

```bash
    doxygen doc/Doxyfile
    firefox doc/html/index.html
```

## Author

StatAnaly is solely developed by Ansel Blumers.

## Copyright

> Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0

> Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

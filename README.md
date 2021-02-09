# UKFilter_STM32F4
## About
Unscented Kalman filter for STM32F4 microcontroller implemented in C.  
Project written in VSCode.  
Debugging with the Renode emulator.

## Install
- Install arm-none-eabi-gcc toolchain and add to path. [LINK](https://developer.arm.com/tools-and-software/open-source-software/developer-tools/gnu-toolchain/gnu-rm/downloads/)  
- Install Renode and add to path. [LINK](https://renode.io/)  
- Clone project `git clone https://github.com/IvanVnucec/UKF_STM32F4`
- Go into root directory `cd UKF_STM32F4`
- Pull Git submodules `git submodule update --init --recursive`
- Open with VSCode `code .`

## Building
- Position yourself in root dir and run `make`

## Testing
Position yourself in root dir and run `make start_renode`


## Acknowledge
Special thanks to [ivo-georgiev](https://github.com/ivo-georgiev) for his [UKF Library](https://github.com/ivo-georgiev/ukfLib).

## License
MIT license. See LICENSE file.

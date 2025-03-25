#!/bin/bash

setup_oneapi() {
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
    |   gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
    | sudo tee /etc/apt/sources.list.d/oneAPI.list

    sudo apt update
    sudo apt install intel-basekit intel-hpckit
}

setup_fortran_project(){
    echo "Setting up fortran_project script..."
    mkdir -p ~/.local/bin
    script="$(curl https://raw.githubusercontent.com/ipqa-research/fortran-setup/main/fortran_project)"
    echo "$script" > ~/.local/bin/fortran_project
    chmod +x ~/.local/bin/fortran_project
}

install_packages() {
    # Update apt repositories
    sudo apt update && sudo apt upgrade

    sudo apt install \
        git \
        python3-pip python3-venv pipx \
        gfortran \
        libblas-dev liblapack-dev \
        gdb \
        fzf

    # Install fortran language server, fprettify and flinter
    packages=( fortls findent flinter ford fpm fypp fprettify )
    for package in ${packages[@]}; do
        pipx install $package --force
    done
}


ask_yn() {
    echo "$1 [y/n]"
    read var
    [ $var = "y" ]
}

setup_all() {
    setup_fortran_project
    setup_oneapi
    install_packages
}

echo "
Software Bootstrap necessary to run Fortran code on Linux

This script will install:
- Mathematical libraries (LAPACK) for solving linear systems
- Compiler and debugger gfortran
- OPTIONAL: Intel OneAPI HPC Toolkit (Intel compilers, weighs about 10GB)
- Development tools:
    - pipx: Python program manager
    - fortls: To display errors and Fortran syntax
    - findent: To format code (order indentations)
    - flinter: To perform code analysis
    - fpm: Code package manager, for generating and running projects
    - Fortran project generator script fortran_project
"

if ask_yn "Install everyhing?"; then
    setup_all
else
    ask_yn "Install oneAPI" && setup_oneapi
    ask_yn "Install tools?" && install_packages
    ask_yn "Install fortran_project" && setup_fortran_project
fi

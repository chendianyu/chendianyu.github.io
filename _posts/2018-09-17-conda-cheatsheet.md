---
title: Conda cheatsheet
description: Conda 简介
categories:
 - Conda
tags:
 - Conda
 - software
---

# Conda
基本操作  
```shell
conda --version
conda -V 
conda info
conda update conda
```

# Environment  
Conda `environment` 指的是包含所安装的特定 conda 包的路径，默认的环境名称为 `base`。当你改变某一环境时，其他的环境并不会受到影响。可以通过 `activate` 和 `deactivate` 实现环境的切换。另外还可以通过将 `environment.yaml` 拷贝给他人实现环境的共享    
  
## create an environment  
```shell
conda create --name <myenv>    # create an environment called myenv
```
By default, environments are installed into the `envs` directory in your conda directory  
  
## remove an environment 
```shell
conda remove --name <myenv> --all
```
  
## activate and deactivate an environment  
```shell  
source activate <myenv>  
source deactivate       # Change your current environment back to the default (base)
```

## list all your environments  
```shell 
conda info --envs
conda env list
```
  
# Packages  
```shell
conda install <package>
conda search <package>    # Check to see if a package you have not installed is available from the Anaconda repository
conda list                # Check to see packages in the environment
conda update <package>    # update a package
conda remove <package>    # remove a package
# 以上操作均在当前env中进行 如果想指定某个环境 则加上 -n <env> 即可
```
  
# Configuration  
## .condarc conda configuration file
`.condarc` 文件可以包含多种设置:  
* conda 从哪里搜索包 (Where conda looks for packages)  
* conda 是否以及如何使用代理服务器 (If and how conda uses a proxy server)  
* 创建新环境时同时包含的默认的包或特征 (Default packages or features to include in new environments)  
  
创建或修改 `.condarc` 文件, 可以使用命令 `conda config` 或通过文本编辑器直接操作 `.condarc` 文件，并将其保存至用户主目录（第一次使用 conda config 时会默认创建） or root directory  
  
## channels  
使用 `defaults` 可以自动包含所用默认的通道  
如果仅希望在某个环境中进行通道设置，可以将 `.condarc` 文件放在该环境的 root  directory 下, e.g. `~/miniconda3/envs/test/.condarc`，即仅作用于 test 环境；或者在使用 `conda config` 时使用 `--env` 参数  
  
参考配置命令：  
```shell
conda config --add channels conda-forge    # 后写入的通道将会放在最前面
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes
```

`.condarc` 文件示范  
```
channels:
  - https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
  - https://mirrors.ustc.edu.cn/anaconda/pkgs/main/
  - https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
  - conda-forge
  - defaults
```
  
## Update conda automatically
When `auto_update_conda: True` (default), conda updates itself any time a user updates or installs a package in the root environment. When `False`, conda updates itself only if the user manually issues a `conda update` command  
  
## Show channel URLs (`show_channel_urls`)  
Show channel URLs when displaying what is going to be downloaded and in `conda list`. The default is `False`  
  
## Add pip as Python dependency (`add_pip_as_python_dependency`)
Add pip, wheel and setuptools as dependencies of Python. This ensures that pip, wheel and setuptools are always installed any time Python is installed. The default is `True`  
  
## Always add packages by default (`create_default_packages`)
When creating new environments, add the specified packages by default. The default packages are installed in every environment you create. You can override this option at the command prompt with the `--no-default-packages` flag. The default is to not include any packages  
```shell
create_default_packages:
  - pip
  - ipython
  - scipy=0.15.0  
```

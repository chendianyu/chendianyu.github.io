---
title: Linux 软链接和硬链接
description: 
categories:
 - Linux
tags:
 - Linux
 - concept
---

理解 Linux 中的软链接和硬链接。  

<!-- more -->

# Linux 文件系统  
Linux 有着极其丰富的文件系统，大体上可分如下几类：  
* 网络文件系统，如 nfs、cifs 等；  
* 磁盘文件系统，如 ext4、ext3 等；  
* 特殊文件系统，如 proc、sysfs、ramfs、tmpfs 等  
  
实现以上这些文件系统并在 Linux 下共存的基础是 `Linux VFS（Virtual File System` 又称 `Virtual Filesystem Switch`，即虚拟文件系统。VFS 作为一个通用的文件系统，抽象了文件系统的四个基本概念：文件、目录项 (dentry)、索引节点 (inode) 及挂载点，其在内核中为用户空间层的文件系统提供了相关的接口（下图）。VFS 实现了 open()、read() 等系统调并使得 cp 等用户空间程序可跨文件系统。VFS 真正实现了在 Linux 中除进程之外一切皆是文件  
![vfs](/img/2018-06-24-soft-and-hard-link/vfs.png)  
  
文件包含文件名和数据，而后者在 Linux 中又分为`用户数据（user data）`和`元数据（metadata）`。用户数据，即文件数据块 (`data block`)，是记录文件真实内容的地方；而元数据则是文件的附加属性，如文件大小、创建时间、所有者等信息。  
在 Linux 中，元数据中的 `inode` 号（inode 是文件元数据的一部分但其并不包含文件名，inode 号即索引节点号）才是文件的唯一标识，而不是文件名。文件名仅是为了方便人们的记忆和使用，系统或程序通过 inode 号寻找正确的文件数据块。下图展示了程序通过文件名获取文件内容的过程  
![inode](/img/2018-06-24-soft-and-hard-link/inode.png)  
  
为解决文件的共享使用，Linux 系统引入了两种链接：硬链接 (`hard link`) 与软链接（又称符号链接，即 `soft link` 或 `symbolic link`）。链接为 Linux 系统解决了文件的共享使用，还带来了隐藏文件路径、增加权限安全及节省存储等好处。  
  
# 硬链接  
硬链接会创建独立的虚拟文件，其中包含了原始文件的信息及位置，但是它们从根本上而言是**同一个文件**，因此它们对应的是**相同的 inode 号**。硬链接相当于使同一个文件具有多个别名，引用硬链接文件等同于引用了源文件  
Usage: `ln <oldfile> <newfile>`  
由于硬链接是有着相同 inode 号仅文件名不同的文件，因此硬链接存在以下几点特性：  
* 文件有相同的 inode 及 data block；  
* 只能对已存在的文件进行创建；  
* 不能跨文件系统进行硬链接的创建；  
* 不能对目录进行创建，只可对文件创建；  
* 删除一个硬链接文件并不影响其他有相同 inode 号的文件  
  
![link](/img/2018-06-24-soft-and-hard-link/link.png)  
  
# 软链接  
软链接与硬链接不同，它是一个事实存在的文件，只不过用户数据块中存放的内容是另一文件的路径名的指向。  
Usage: `ln -s <oldfile> <newfile>`  
软链接的创建与使用没有类似硬链接的诸多限制：  
* 软链接有自己的文件属性及权限等；  
* 可对不存在的文件或目录创建软链接；  
* 软链接可交叉文件系统；  
* 软链接可对文件或目录创建；  
* 创建软链接时，链接计数 i_nlink 不会增加；  
* 删除软链接并不影响被指向的文件，但若被指向的原文件被删除，则相关软连接被称为死链接（即 `dangling link`，若被指向路径文件被重新创建，死链接可恢复为正常的软链接）  
  
# REF  
1. https://www.ibm.com/developerworks/cn/linux/l-cn-hardandsymb-links/index.html  
2. Linux 命令行与 shell 脚本编程大全（第三版）

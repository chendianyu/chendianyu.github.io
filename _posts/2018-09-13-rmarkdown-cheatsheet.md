---
title: R Markdown cheatsheet
description: R Markdown 简介
categories:
 - R
tags:
 - R
 - Markdown
---

# Introduction
首先来理清几个基本的概念。`markdown` 是一种轻量级的标记语言，适合记录技术文档以及笔记等；`Pandoc` 包用于不同格式标记语言文档之间的转换，例如将 markdown 文档转换成 pdf，HTML 等；`knitr` 包能将代码及运行结果等嵌入到 markdown 中，并将 R Markdown 转换成 markdown；`rmarkdown` 包则是站在 knitr 和 Pandoc 的基础上，完成了整个任务。  
  
# Basics
R Markdown 文档包含3个基本的组成部分：`metadata`（可选）, `text`, and `code`  
metadata 部分用一对 `---` 包围，语法格式为 YAML，所以也被称为 YAML 元数据或 YAML 前辅文。  
```rmarkdown
---
title: "Hello R Markdown"
author: "Awesome Me"
date: "2018-02-14"
output: html_document
---
```  

text 部分的语法格式为 Pandoc's Markdown  
几个注意的点包括：  
* 下标用一对 `~`，上标则用一对 `^`  
* 行内代码，如果代码里存在反引号，有n个，那么外面包裹它的反引号则需要n+1个，且与文本之间有空格  
* 脚注放在 `^[]` 内  
* 数学表达式参照 LaTeX 语法  
  
code 部分有两种形式：  
* 代码块，用3个反引号 ```` ```{r} ```` 起始，其中 r 表示使用的编码语言的名称，另外还可以在花括号中设置代码块选项，如字体大小等，最后用3个反引号结束  
* 行内代码，则以 `` `r `` 起始，以 `` ` `` 结束，直接运行代码并在最终文档中展示结果  
  
rmarkdown 包中提供了两种输出格式：documents 和 presentations。当在 YAML metadata 中指定输出格式时，如果是调用的第三方包的，则还需要同时注明是哪个包。 e.g.: `output: tufte::tufte_html`。
  
# Code Chunks
在 Rstudio 中，可以手动输入 ```` ```{r} ```` 和 ```` ``` ```` 表示代码块，也可以用快捷键 `Cmd/Ctrl-Alt-I`  
可以跟运行其他 r 代码一样，通过快捷键 `Cmd/Ctrl-Enter` 运行 chunk 中的代码，也可以用 chunk 专属的快捷键 `Cmd/Ctrl-Shift-Enter` 运行 chunk 中所有的代码  
chunk header 包含多种参数可供定义，接下来让我们看看它们的用法  
  
## chunk name
可以在 header 中设置代码块的名称 ```` ```{r <name>} ````，这样便于导航，以及有代码产生的图方便在其他地方继续使用  
`setup` 是一个特殊的名字，当中的代码将会在其他所有代码执行之前便自动运行一次  
  
## chunk options  
`knitr` 提供了约60种选项用于自定义我们的代码块，具体的可以参考 <http://yihui.name/knitr/options/>  
最主要的几种是对代码运行以及结果展示的控制：  
* `eval = FALSE`，**阻止代码的运行**，常用于代码的展示  
* `include = FALSE`，**运行代码，但不在最终文档中展示代码和结果**，用于一些设置  
* `echo = FALSE`，**运行代码，不显示代码，但显示结果**，适用于不需要知道背后运行的具体代码的情况  
* `message = FALSE` or `warning = FALSE`，**不在最终结果中显示 message 和 warnings**  
* `results = 'hide'` and `fig.show = 'hide'`，**分别隐藏运行结果和图片**  
* `error = TRUE`，**即便代码报错，仍进行渲染**，适合调试。默认的 `error = FALSE`，只要存在错误，就会导致 knit 中止  

下图是上述几种选项的总结  
![code_chunk_option](/img/2018-09-13-rmarkdown-cheatsheet/code_chunk_option.png)  
  
## Cache  
每个文档的编织过程都是全新的，但有时候有些计算花费时间过长，老是重新跑也是件痛苦的事，可以通过 `cache = TRUE` 将代码的输出存储到硬盘上。在之后运行时，knitr 将会检查代码是否发生变化，如果没有就会直接使用该结果。需要注意的是缓存**仅基于该块代码本身，不包括其依赖关系**，所以即便其依赖的代码发生了变化，它也不会重跑。为了解决这一问题，可以再加上 `dependson` 参数，包含每个所依赖代码块名字的**字符向量**， e.g. `dependson = "raw_data"`，从而确保当依赖的代码发生变化时，也能够察觉  

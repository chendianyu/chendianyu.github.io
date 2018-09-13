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
首先来理清几个基本的概念。`markdown` 是一种轻量级的标记语言，适合记录技术文档以及笔记等；`Pandoc` 包用于不同格式标记语言文档之间的转换，例如将 markdown 文档转换成 pdf，HTML 等。`knitr` 包能将代码等嵌入到 markdown 中，并将 R Markdown 转换成 markdown；`rmarkdown` 包则是站在 knitr 和 Pandoc 的基础上，完成了整个任务。  
  
# Basics
R Markdown 文档包含3个基本的组成部分：`metadata`, `text`, and `code`  
metadata 部分用一对 `---` 包围，语法格式为 YAML，所以也被称为 YAML 元数据或 YAML 前辅文。  
```r markdown
---
title: "Hello R Markdown"
author: "Awesome Me"
date: "2018-02-14"
output: html_document
---
```  
text 部分的语法格式为 Markdown，其中计算机代码部分有两种格式：  
* 代码块用三个反引号 `\`\`\`{r}` 

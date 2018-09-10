---
title: Pandas 使用手册
description: Pandas 作为 Python 中的第三方数据处理包，给我们的数据分析带来了很大的便利，在此将当中常用到的数据处理函数进行简单的整理，方便日后的查阅。
categories:
 - RNA-seq
tags:
 - RNA-seq
 - bioconductor
---

# Data Loading  
* `pd.read_csv()/read_table()`  
  * `sep`： 默认分别为 `','` and `'\t'`，也可以是正则表达式  
  * `header`： 默认第一行为 header；`None` 则无 header，显示时为 range(n)； `<int>` 表明第 int+1 行做 header；如果是 `list of int`，则生成multi-index  
  * `names`: array-like (List of column names to use)，自定义列名  
  * `index_col`: int or name of column(s), 指定作为index的列；多列就形成 multi-index  
  * `skiprows`: list-like or integer, 跳过指定行  
  * `comment`: 注释行，从而忽略 e.g. `comment='#'`  
  * `na_values`: scalar, str, list-like, or dict, 认为是 nan 的值  
  
# Data Saving  
* `df.to_csv()`  
  * `path`  
  * `sep`  
  * `index/header=True`: 是否输出 index 和 columns 名称  
  * `na_rep`: string, 缺失值的输出结果, default to `''`  
  * `columns`: list of columns to write  

# Basic attributes  
* shape (返回一个 tuple)
* dtypes  
* index  
* columns  
* values  
* index/columns.name (针对普通的， type是 `str`)  
* index/columns.names (针对MultiIndex，普通的也可以用，只是type为 `pandas.core.indexes.frozen.FrozenList`)  
* index/columns.is_unique: 返回唯一值  
* iloc/loc: slicing, indexing  
(`pd.Index()`, `pd.Series()` 等有 `name` 参数，`pd.MultiIndex()` 有 `names` 参数，`name` 参数也有，但只是为了兼容)  
  
# Data viewing
```python
df.head(<num>)    # default to 5
df.tail(<num>)
```

# Index reconstruction  
* `reindex()`: 改变index的组成以及顺序等  
  * `index/columns`: list of labels, 按给出的顺序以及名字重构  
  * `method`: `{'ffill',...}`, 指定填充引入的缺失值的方法  
  
* `set_index(<column name or list of column names>)`: 选取列作为新的index  
  * `drop=True`: 是否在values中仍保留这些列  
* `reset_index()`: 默认将所有的index level转回成column，处于最前列 (更详细的用法见官方文档)  
  
* `sort_index()`: 对index或columns进行排序  
  * `axis=0`  
  * `level`: multi-index中使用，对level排序  
  * `ascending=True`  
  
* `drop()`: 舍弃指定index或columns  
  * `labels`: single label or list-like  
  * `axis`  
  * `level`  
  * `inplace=False`  
  
* `rename()`: 对index或columns进行重命名，也可以以`map`方法作用于index或columns  
  * `index/columns`: 可以是函数，或者字典，其中key为原名，value为新名  
  * `inplace=False`  
  
# Value reconstruction  
* `sort_values()`  
  * `by`: str or list of str, 指定用于排序的key  
  * `axis`  
  * `ascending=True`  
  * `na_position`: ['first', 'last'], 默认缺失值排最后，即last  
  * `inplace=False`  
  
* `duplicated()`: Return boolean Series denoting duplicate rows  
  * `subset`: column label or sequence of labels, 根据指定列来判断是否重复  
  * `keep`:  
    * `first`: Mark duplicates as True except for the first occurrence  
    * `last`: Mark duplicates as True except for the last occurrence  
    * `False`: Mark all duplicates as True  
* `drop_duplicates()`: 删除重复行，等价于 `df[~df.duplicated()]`  
  * `subset`  
  * `keep`  
  * `inplace`  
  
* `replace()`: 替换值, 可以是 原值（多个值就用list），新值； 如果是不同值分别替换，每一对可以是list或dict  
```python  
data.replace(-999, np.nan)
data.replace([-999, -1000], np.nan)
data.replace([-999, -1000], [np.nan, 0])
data.replace({-999: np.nan, -1000: 0})
```  
  
# Missing value  
* `dropna()`: 根据缺失值情况删除行或列  
  * `axis`  
  * `how`: {'all', 'any'}, 默认为any，即存在缺失值就删除整行或整列  
  * `thresh`: int, 指定label中有效值少于int，则删除该label  
  * `inplace`  
  
* `fillna()`  
  * `value`: scalar（所有缺失值均替换成该值）或 dict（不同的列对应不同的替换值，貌似只有对axis=0可行？）  
  * `method`: {'ffill',...}, 填补缺失值的方法  
  * `limit`: int, If method is specified, this is the maximum number of consecutive NaN values to forward/backward fill 
  * `inplace`  
  
# Merge, join and concatenate  
* `pd.merge()` (一般原来的index会被抛弃，除了`left_index`和`right_index`均为`True`时)  
  * `left/right`: DataFrame  
  * `how`: {'inner', 'outer', 'left', 'right'}, 默认inner  
  * `on`: label or list, 指定作为合并的key的label  
  * `left_on/right_on` (选作key的列均会在values中出现，所以可能需要用 `drop()` 剔除)  
  * `left_index/right_index=False`  
  * `suffixes`: 2-length sequence, 默认('\_x', '\_y')  
  * `sort=False`: 是否根据合并时使用的key来进行字典排序  
(df.merge()类似，df.join()默认根据index合并，且默认left join)  
  
* `pd.concat()`  
  * `obj`: a sequence of Series, DataFrame  
  * `axis`  
  * `join`: {'outer', 'inner'}, 默认outer  
  * `join_axes`: list of labels, 结果中保留的labels, 是与axis指定的不同axis上的label  
  * `keys`: sequence, 长度与`obj`一致，从而构造出multi-index，并作为最外层的level (也算是一种避免原objs的index的重复导致结果混乱的方法)  
  * `names`: list, 生成的multi-index各level的name  
  * `ignore_index=False`: 是否忽略原来的index，而使用0,...,n-1  
  * `verify_integrity=False`: 检查结果中index是否存在重复，若重复则raise exception  
  
# Reshaping  
* `stack()`: Pivot a level of the (possibly hierarchical) `column labels`, returning a DataFrame (or Series in the case of an object with a single level of column labels) having a hierarchical index with a new `inner-most` level of row labels. The level involved will `automatically get sorted`  
  * `level`: int, string, or list of these, default last level, stack使用的level  
  * `dropna=True`: Whether to drop rows in the resulting Frame/Series with no valid values  
* `unstack()`: Pivot a level of the (necessarily hierarchical) `index labels`, returning a DataFrame having a new level of column labels whose `inner-most` level consists of the pivoted index labels. If the index is not a MultiIndex, the output will be a Series (the analogue of stack when the columns are not a MultiIndex). The level involved will `automatically get sorted`  
  * `level`  
  * `fill_value`: replace NaN with this value if the unstack produces missing values  
  
* `pivot()`: Reshape data based on `column values`. Uses unique values from index / columns to form axes of the resulting DataFrame  
  * `index=None` : string or object, column name to use to make new frame's index. Default to uses existing index， 且**不能有重复**(与`pivot_table()`的区别)  
  * `columns=None` : string or object, column name to use to make new frame's columns  
  * `values=None` : string or object, Column name to use for populating new frame's values. By default, all remaining columns will be used and the result will have hierarchically indexed columns, 除了被选为index和columns的列外，其余列成为多维列最外层的level  
```python
>>> df = pd.DataFrame({'foo': ['one','one','one','two','two','two'],
                       'bar': ['A', 'B', 'C', 'A', 'B', 'C'],
                       'baz': [1, 2, 3, 4, 5, 6],
                       'bac': [-1,-2,-3,-4,-5,-6]})
>>> df.pivot(index='foo', columns='bar', values='baz') # 相当于df.pivot(index='foo', columns='bar')['baz']
     A   B   C
one  1   2   3
two  4   5   6
>>> df.pivot(index='foo', columns='bar')
    bac         baz
bar A   B   C   A   B   C
foo                     
one -1  -2  -3  1   2   3
two -4  -5  -6  4   5   6  
```  
  
* `pd.melt()`: merges multiple columns into one, producing a DataFrame that is longer than the input  
  * `frame`: DataFrame  
  * `id_vars` : tuple, list, or ndarray, column(s) to use as identifier variables, 用来做group indicator的列，不会被合并； 即便只选了一列，也得是list形式的； 也可以什么都不选  
  * `value_vars` : tuple, list, or ndarray, column(s) to unpivot, 即要合并的列. Default to uses all columns that are not set as `id_vars`， 没选的列就被扔掉了
  
# Groupby  
`groupby(by=None, axis=0, level=None, as_index=True, sort=True, group_keys=True, squeeze=False, **kwargs)`  
* `by`: mapping, function, str, or iterable  
  * A str or list of strs may be passed to group by the columns  
  * If a function, it's called on each value of the object's index (把具有相同函数返回结果的作为一组)  
  * If a `dict` or `Series` is passed, the Series or dict VALUES will be used to determine the groups (the Series' values are first aligned) (index对应相同VALUE的就归为一组，例如`{'A': 'vowel', 'B': 'consonant', 'C': 'consonant'}`, B和C的就作为一组)  
  * If an `ndarray` is passed, the values are used as-is determine the groups  
* `as_index` : boolean, default True. For aggregated output, return object with group labels as the index；否则就把index转成values最前面的列  
* `sort` : boolean, default True. Sort group keys. Note this does not influence the order of observations within each group. groupby preserves the order of rows within each group  

* `grouped.aggregate()` 或 `grouped.agg()`: 可以针对每一列实现更加灵活地使用函数  
  * 单个函数或string of function name  
  * list of functions，从而对每列数据实现多个统计  
  * dict of column name:functions (or list of functions)，不同列可以给定不同的函数，结果中也只会包含这些列，其他的省略了；不可以直接多列放在一起，因为会判定成keyerror  
  
* `grouped.filter(func, dropna=True, *args, **kwargs)`: Return a copy of a DataFrame excluding elements from groups that do not satisfy the boolean criterion specified by func. 即利用func作用于各分组，得出结果为True的分组，然后利用该值去对原DataFrame进行slicing  
  * `func`: Function to apply to each subframe. Should return True or False  
  * `dropna`: Drop groups that do not pass the filter; if False, groups that evaluate False are filled with `NaN`  
```python
>>> rng = np.random.RandomState(0)
>>> df = pd.DataFrame({'key': ['A', 'B', 'C', 'A', 'B', 'C'],
                   'data1': range(6),
                   'data2': rng.randint(0, 10, 6)},
                   columns = ['key', 'data1', 'data2'])
>>> def filter_func(x):
        return x['data2'].std() > 4
>>> print(df.groupby('key').filter(filter_func))
  key data1 data2
1 B   1     0
2 C   2     3
4 B   4     7
5 C   5     9
```  
  
* `grouped.transform(func, *args, **kwargs)`: Call function producing a like-indexed DataFrame on each group and
return a DataFrame having the same indexes as the original object filled with the transformed values  
```python
>>> df.groupby('key').transform(lambda x: x - x.mean())
  data1 data2
0 -1.5  1.0
1 -1.5  -3.5
2 -1.5  -3.0
3 1.5   -1.0
4 1.5   3.5
5 1.5   3.0
``` 
  
* `grouped.apply(func, *args, **kwargs)`: apply an arbitrary function to the group results. The function should take a `DataFrame`, and return either a `Pandas object` (e.g., DataFrame, Series) or a `scalar`; the combine operation will be tailored to the type of output returned  
  
# Selection
## selection by label .loc
The `.loc` attribute accept the following valid inputs:  
  * A single label, e.g. `5` or `'a'`  
  * A list or array of labels `['a', 'b', 'c']`  
  * A slice object with labels `'a':'f'` (**both the start and the stop are included**)  
  * A boolean array  
  * A callable  
```python
# 注意label的匹配，尤其是通过boolean array进行索引时
>>> df1 = pd.DataFrame(np.random.randn(6,4),
                       index=list('abcdef'),
                       columns=list('ABCD'))
>>> df1
>>> 
    A         B         C         D
a  0.132003 -0.827317 -0.076467 -1.187678
b  1.130127 -1.436737 -1.413681  1.607920
c  1.024180  0.569605  0.875906 -2.211372
d  0.974466 -2.006747 -0.410001 -0.078638
e  0.545952 -1.219217 -1.226825  0.769804
f -1.281247 -0.727707 -0.121306 -0.097883
>>> df1.loc['a'] > 0
>>> 
A     True
B    False
C    False
D    False
Name: a, dtype: bool
>>> df1.loc[:, df1.loc['a'] > 0]
>>> 
     A
a  0.132003
b  1.130127
c  1.024180
d  0.974466
e  0.545952
f -1.281247
>>> df1.loc[:, df1.loc['a'] > 0]
>>> IndexingError: Unalignable boolean Series provided as indexer (index of the boolean Series and of the indexed object do not match
```
  
## selection by position .iloc
The `.iloc` attribute accept the following valid inputs:  
  * An integer e.g. `5`  
  * A list or array of integers `[4, 3, 0]`  
  * A slice object with ints `1:7`  
  * A boolean array  
  * A callable  
```python
# 由于boolean array通常是带有index的，所以一般不会直接做iloc的输入
# 如果非要做输入，就得使用其values属性
>>> df2 = pd.DataFrame(np.random.randn(4,4),
                      index=list('abcd'),
                       columns=list('ABCD'))
>>> df2
>>>   A   B   C   D
a   1.561034  1.004562  -0.404225   -1.232107
b   -0.804208   -1.746975   0.643511  0.156833
c   -0.827571   0.480988  -0.977209   1.430038
d   -0.395301   0.974563  1.375785  -1.065290 
>>> df2.iloc[df2['A']>0]
>>> ValueError: iLocation based boolean indexing cannot use an indexable as a mask
>>> 
>>> df2.iloc[(df2['A']>0).values]
>>> A   B   C   D
a   1.561034  1.004562  -0.404225   -1.232107
>>> df2.iloc[:, (df2['A']>0).values]
>>>   A
a   1.561034
b   -0.804208
c   -0.827571
d   -0.395301
```

## Boolean indexing  
The operators are:  
  * `|` for or  
  * `&` for and  
  * `~` for not  
(表达式都用括号括起来以防止歧义，e.g. `s[~(s>0)]`)  



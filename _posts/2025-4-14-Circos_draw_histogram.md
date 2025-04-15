---
layout: post
title: 使用circos画变异分布直方图
date: 2025-04-14 22:00
tags: [SV,VISUAL]
toc: true
---

## Installation
---
Circos主要使用Perl编写，依赖了大量的Perl包，官网有详细的安装教程，本文不再介绍，贴上部分链接。

Perl & Module Install :  [perl and perl modules](https://circos.ca/documentation/tutorials/configuration/perl_and_modules/)




Circos install: [Install circos](https://circos.ca/documentation/tutorials/configuration/installation/)



## Description
---
Circos的圈图基本组成：最外圈染色体（本人这样干，可以不这样）称为：karyotype，里面的一圈为一种数据，称为一个track。

Circos作图完全依赖配置文件，在配置文件中指定需要可视化的数据，可视化的类型，以及图像的各种参数，见下部分。



## Configuration
---
基础配置文件描述

```
karyotype = "染色体文件"

<ideogram>
各种基础配置
</ideogram>

<image>
图片配置,包括图片大小,像素,保存位置等
</image>

<plots>
所有图像
<plot>
单个图像
</plot>
</plots>

<links>
线条
</links>
```

## Data preperation
---


### 准备染色体文件

染色体文件即上述的karyotype.txt , 该文件基本格式如下

```
chr - CHRNAME CHRLABEL START END COLOR 
column 1: 表示该行记录是染色体
column 2: - 必要的
column 3: 染色体名字,看你自己
column 4: 染色体标签,可以理解为最后图上显示的染色体名字
column 5: 染色体起始位置
column 6: 染色体终点
column 7: 颜色,我这里指定chr1,chr2 etc是因为有内置的颜色名为chrN,后面可以自己指定

chr - chr1 1 0 5000000 chr1
chr - chr2 2 0 10000000 chr2 
chr - chr3 3 0 20000000 chr3 
chr - chr4 4 0 50000000 chr4 
chr - chr5 5 0 100000000 chr5
```

一般来说，确定染色体后，就不需要再修改了，这个可以一直使用，本文只介绍了单个圈图的情况，还有其他多图情况可以参看官方文档。



### 准备数据文件

前面提到，内圈的数据都是一个track对应一个数据文件，那这个文件格式如何，怎么准备呢?

本文主要描述变异Histogram怎么做，其他类型请自行查阅文档。


- 获取基因组范围

我们需要准备一个基因组范围文件，以制表符分隔，格式为: CHR \t START(0) \t END 

如果你使用RefSeq参考基因组，你可以下载其sequence_report.tsv 使用awk工具提取，最终结果应如下

```

染色体名字 \t 0 \t 染色体终点  
CHR1 \t 0 \t 5000000

```


- 获取变异位点

如果是snp、indel或者SV(单纯看分布，不需要SVTYPE)

```shell
bcftools query -f '%CHROM\t%POS' your_vcf_file > variant_site.txt
```

如果是SV并且需要按照SVTYPE分组，作堆叠直方图

```shell
bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE' your_vcf_file > variant_site.txt
```

3. 获取分箱文件

直方图，需要一个范围并且统计范围内变异的数量，这个过程称为分箱，在准备好了上述两文件后，使用下列Python脚本获取分箱文件

```
import sys
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("genome", help="genome file: chr, start, end")
parser.add_argument("variant", help="variant file: chr, pos[, type]")
parser.add_argument("bin_size", type=int)
parser.add_argument("--typed", action="store_true", help="3-column input with types")
args = parser.parse_args()

# Step 1: 读取基因组范围
genome_df = pd.read_csv(args.genome, sep='\t', header=None, names=['chr', 'start', 'end'])

# Step 2: 读取变异信息
if args.typed:
    var_df = pd.read_csv(args.variant, sep='\t', header=None, names=['chr', 'pos', 'type'])
else:
    var_df = pd.read_csv(args.variant, sep='\t', header=None, names=['chr', 'pos'])
    var_df['type'] = 'COUNT'  

result_list = []

for _, row in genome_df.iterrows():
    chr_name, chr_start, chr_end = row['chr'], row['start'], row['end']
    # 取出位于该染色体的所有变异
    chr_data = var_df[var_df['chr'] == chr_name].copy()
    if chr_data.empty:
        continue

    # 获取该染色体上每个bin的范围
    bins = np.arange(chr_start, chr_end , args.bin_size)
    if bins[-1] < chr_end:
        bins = np.append(bins,chr_end + 1)
   
    chr_data['bin'] = pd.cut(chr_data['pos'], bins=bins,  right=False,include_lowest=True)
    # 统计每个 bin 的每种 type 数量
    bin_counts = chr_data.groupby(['bin', 'type'],observed=False).size().unstack(fill_value=0)
    # 按 bin 拆分出 bin_start / bin_end
    bin_counts = bin_counts.reset_index()
    bin_counts['bin_start'] = bin_counts['bin'].map(lambda x : int(x.left))
    bin_counts['bin_end'] = bin_counts['bin'].map(lambda x : int(x.right) - 1)
    
    bin_counts.insert(0, 'chr', chr_name)
    bin_counts.drop(columns=['bin'], inplace=True)

    result_list.append(bin_counts)

# 合并所有染色体
all_bins = pd.concat(result_list, ignore_index=True)

# 排序列（按平均值降序）
if args.typed:
    type_cols = all_bins.columns.difference(['chr', 'bin_start', 'bin_end'])
    mean_order = all_bins[type_cols].mean().sort_values(ascending=False).index
    all_bins = all_bins[['chr', 'bin_start', 'bin_end'] + list(mean_order)]
    type_col_name = ','.join(mean_order)
    all_bins[type_col_name] = all_bins[mean_order].astype(str).agg(','.join, axis=1)
    all_bins = all_bins[['chr', 'bin_start', 'bin_end',type_col_name]]
    
else:
    all_bins = all_bins[['chr', 'bin_start', 'bin_end', 'COUNT']]

# 输出
all_bins.to_csv(sys.stdout, sep='\t', index=False)
```


保存上述脚本文件并运行

```shell
python 上述python文件路径 基因组范围文件 变异位点文件 分箱大小 --typed > 分箱文件
```

注意，针对结构变异的堆叠直方图，本脚本在输出分箱文件时，第一行是表头，方便用户查看各列对应的变异种类，作图之前需要将其删除。
## 作图
---
有了染色体文件和分箱文件后，现在完善配置文件，假设你项目文件结构如下

```
project
├── data
│   ├── bin_data.txt
│   ├── karyotype.txt
└── circos.conf
```


circos.conf就是最重要的配置文件，其内容应该如下

```txt
katyotype = ./data/karyotype.txt
# 默认显示所有染色体,如果有其他设置需要将这个设置成NO
chromosomes_display_default = yes 

<ideogram>
	# 表示每条染色体间的距离,r表示相对距离
	<spacing>
		default = 0.005r 
	</spacing>
	
	radius = 0.9r
	thickness = 50p
	fill = yes
	stroke_thickness = 2p
	stroke_color = dgrey
	
	show_label = yes
	lab_font = deafult
	label_radius = 1r + 75p
	label_size = 30
	lable_parallel = yes
</ideogram>

<image>
	# 图片配置,没有自己配置的情况下,使用绝对路径导入位于 circos安装目录/etc/下的image.conf
	<<include circos安装目录/etc/image.conf>>
</image>

<plots>
	<plot>
		type = histogram
		file = ./data/bin_data.txt
		thickness = 0.5p
		color = black
		# max代表你想在图中显示的最大值,根据你每个窗口的值来看
		max = 20000
		# r0表示该track的内半径,单位r表示相对于这个圆的大小
		r0 = 0.85r 
		# r1表示该track的外半径,单位r表示相对于这个圆的大小
		r1 = 0.95r
		# 注意,我这里使用的SV的数据,并且打算做堆叠柱状图,我的堆叠中有4个种类,所以这里有四种颜色,如果只有一种的只需要一种颜色,这里的颜色的定义在Circos安装目录下etc下的colors.conf中,需要导入
		fill_color = vvdblue,dblue,lblue,vvblue
	</plot>
</plots>


<<include circos安装目录/etc/colors_fonts_patterns.conf>>
<<include circos安装目录/etc/colors.conf>>
<<include circos安装目录/etc/housekeeping.conf>>

```


准备好上述工作后，只需要在项目根目录下运行如下命令，将会得到circos.png 和 circos.svg两种格式图片
```shell
circos --conf circos.conf
```




![circos.svg](https://raw.githubusercontent.com/kunmonster/note_pic/main/note/202504142100655.svg)


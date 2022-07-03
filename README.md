# Snakemake workflow: rna-seq-star/hisat2/salmon-deseq2

- 使用`STAR/hisat2`进行dna比对，并使用`featureCounts`进行计数
- 使用`salmon`进行cdna比对并计数
- 使用`multiqc`生成报告
- 使用`DESeq2`进行差异分析


## Usage

```
# 建立输出文件夹
mkdir /path/to/output

# 生成初始配置文件
python ./setup.py --init

# 创建环境
python ./setup.py -s --conda-create-envs-only

# 启动
python ./setup.py -s
```

修改自<https://github.com/BIMSBbioinfo/pigx_rnaseq>

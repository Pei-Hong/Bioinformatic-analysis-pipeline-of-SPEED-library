#!/bin/bash

file=$1

zcat $file | paste - - - - | grep -v NOT_FOUND | awk -F "\t" '{print $1"\n"$2"\n"$3"\n"$4}' | gzip > ${file%.barcode.fq.gz}.barcode_full.fq.gz

zcat $file | paste - - - - | grep NOT_FOUND | awk -F "\t" '{print $1"\n"$2"\n"$3"\n"$4}' | gzip > ${file%.barcode.fq.gz}.barcode_short.fq.gz # 5hmc_1P.barcode_short.fq.gz

n_full=`zcat ${file%.barcode.fq.gz}.barcode_full.fq.gz | paste - - - - | wc -l`
n_short=`zcat ${file%.barcode.fq.gz}.barcode_short.fq.gz | paste - - - - | wc -l`

echo "n_full:"$'\t'$n_full$'\n'"n_short:"$'\t'$n_short > ${file%.barcode.fq.gz}.barcodeCheck.log


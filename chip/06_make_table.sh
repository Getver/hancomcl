FOLDER='pharmacoFocus/0221'
KEY='0221'

## vcf를 chr 버전으로 바꾼 것(이상한 chr는 딱히 빼지 않음(alt 말한것))
mkdir /disk0/sm/PMDA/${FOLDER}/minjin
zcat /disk0/sm/PMDA/${FOLDER}/09_Export/${KEY}.vcf.gz | grep '^#' > TMP_header.txt
zcat /disk0/sm/PMDA/${FOLDER}/09_Export/${KEY}.vcf.gz | grep -v '^#' | awk '{OF=OFS="\t"; $1="chr"$1; print $0}' | cat TMP_header.txt - > /disk0/sm/PMDA/${FOLDER}/minjin/${KEY}_chr.vcf


# 돌렸더니 뭐 이상한 chr때문에 에러떴는데 그래서 grep 해서 chr이상한건 다 손수 뺐음......(chr27이상 지움, 폴더 자동생성됨, 오래걸림)
java -jar /disk0/sm/PMDA/tool/PharmCAT-2.11.0/pharmcat-2.12.0-all.jar -vcf /disk0/sm/PMDA/${FOLDER}/minjin/${KEY}_chr.vcf --output-dir /disk0/sm/PMDA/${FOLDER}/minjin/pharmCAToutput


## Plink frequency(ped, map 파일 만들었음)

apt-format-result \
--batch-folder /disk0/sm/PMDA/${FOLDER}/06_Genotype \
--cn-region-calls-file /disk0/sm/PMDA/${FOLDER}/05_CN/AxiomCNVMix.cnregioncalls.txt \
--snp-list-file /disk0/sm/PMDA/${FOLDER}/07_SNPolisher/Recommended.ps \
--annotation-file /disk0/sm/chip/library/Axiom_PangenomiX.r1/Axiom_PangenomiX.na36.r1.a1.annot.db \
--snp-identifier-column custom_rsid \
--export-chr-shortname true \
--export-plink-file /disk0/sm/PMDA/${FOLDER}/minjin/${KEY} \
--log-file  /disk0/sm/PMDA/${FOLDER}/minjin/freq.log


mkdir /disk0/sm/PMDA/${FOLDER}/minjin/plink


# bfile 만들고 .frq파일 만들음
plink1.9 --vcf /disk0/sm/PMDA/${FOLDER}/minjin/${KEY}_chr.vcf --allow-extra-chr --make-bed --vcf-half-call m --double-id -out /disk0/sm/PMDA/${FOLDER}/minjin/plink/S1 
plink1.9 --bfile /disk0/sm/PMDA/${FOLDER}/minjin/plink/S1  --allow-extra-chr --freq --make-bed --out /disk0/sm/PMDA/${FOLDER}/minjin/plink/S2





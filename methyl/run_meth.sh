docker run --rm -v /disk2/kb/sm/methyl/meth/01_111:/home/data meth:1.1.0 /home/data/sample_sheet.csv EPICv2
docker run --rm -v /disk2/kb/sm/methyl/meth/01_1:/home/data meth:1.1.0 /home/data/sample_sheet.csv EPICv2

docker run --rm -v /disk2/kb/sm/methyl/meth/01_11:/home/data meth:1.1.0 /home/data/sample_sheet.csv EPICv2
docker run --rm -v /disk2/kb/sm/methyl/meth/01_12:/home/data meth:1.1.0 /home/data/sample_sheet.csv EPICv2
docker run --rm -v /disk2/kb/sm/methyl/meth/01_13:/home/data meth:1.1.0 /home/data/sample_sheet.csv EPICv2

## 결과 알아서 나오게
docker run --rm -v /disk2/kb/Methylation_88ea_DATA/Methylation:/home/data meth:1.1.0 /home/data/sample_sheet.csv EPICv2
docker run --rm -v /disk2/kb/sm/methyl/meth/02:/home/data meth:1.1.0 /home/data/sample_sheet.csv EPICv2

## 내가 들어가서 보게
docker run --rm -it -v /disk2/kb/Methylation_88ea_DATA/Methylation:/home/data meth:test_0.1
docker run --rm -it -v ${PWD}:/home/data sm:meth_0.2


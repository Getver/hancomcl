
LIBRARY=$1
KEY=$2


if [ ${LIBRARY} = "KORV1" ]; then
    bash /code_SM/KORV1.sh ${KEY} ${LIBRARY}
elif [ ${LIBRARY} = "KORV2" ]; then
    bash /code_SM/KORV2.sh ${KEY} ${LIBRARY}
elif [ ${LIBRARY} = "HCLV" ]; then
    bash /code_SM/HCLV.sh ${KEY} ${LIBRARY}
elif [ ${LIBRARY} = "PMDA" ]; then
    bash /code_SM/PMDA.sh ${KEY} ${LIBRARY}
    Rscript /code_SM/plot.R ${KEY} ${LIBRARY}
elif [ ${LIBRARY} = "PangenomiX" ]; then
    bash /code_SM/PangenomiX.sh ${KEY} ${LIBRARY}
    Rscript /code_SM/plot.R ${KEY} ${LIBRARY}
elif [ ${LIBRARY} = "pharmacoFocus" ]; then
    bash /code_SM/pharmacoFocus.sh ${KEY} ${LIBRARY}
    Rscript /code_SM/plot.R ${KEY} ${LIBRARY}
elif [ ${LIBRARY} = "KBA_A" ]; then
    bash /code_SM/KBA_B.sh ${KEY} ${LIBRARY}
elif [ ${LIBRARY} = "KBA_B" ]; then
    bash /code_SM/KBA_B.sh ${KEY} ${LIBRARY}
fi




## test run
# docker run --rm -v ${PWD}:/workdir -v /disk2/hancomcl/KORV1:/disk2/hancomcl/KORV1 sm:chip_0.2 KORV1 TEST
# docker run --rm -v ${PWD}:/workdir -v /disk2/sm/chip/TestData:/disk2/sm/chip/TestData sm:chip_0.2 KORV2 TEST
# docker run --rm -v ${PWD}:/workdir -v /disk2/sm/PMDA/MACROGEN_576/CEL_PangenomiX_96:/disk2/sm/PMDA/MACROGEN_576/CEL_PangenomiX_96 sm:chip_0.2 PangenomiX TEST

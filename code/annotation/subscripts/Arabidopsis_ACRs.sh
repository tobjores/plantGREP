#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### get Arabidopsis DHS union peaks from https://doi.org/10.1093/molbev/msx326 ###
if [ ! -f "${OUTDIR}/Arabidopsis_peaks_union.bed" ]; then
  curl "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mbe/35/4/10.1093_molbev_msx326/8/msx326_supp.zip?Expires=1648882015&Signature=WGDWAUvZUuZqxG1yQdQ230Lc548Kxcp6V2bmTn1lfvBJs2WgRLN-S4ZOrEraalGqtN-FCJ~qSCyO5GJzFuBC0IeecBTePdcoXFJGPB-n-2WaHFybXwqT3t1E8oZEc35~5Y~88U5e4YuD84LO1LHUtd~bdBcy14npTeH~quheFeJ~X52gaCTKZnHy7AJM8bCfTQizzqN0CZuB3yrTos1b8MBbEQctt~4D5Hrcmxp5nosw-EeGo9XPoOduJH1lfs~hVTvsxZ3bsguHWftoxqyw3ZkH2DG1CHL2so5QvZS8-FglnBePdeUMWUbm3P0ebeYF4cD3utW-UYP4JTJoEOlViQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA" \
    > ${TMPDIR}/Alexandre2018_SuppData.zip
  mkdir ${TMPDIR}/Alexandre2018_SuppData
  unzip ${TMPDIR}/Alexandre2018_SuppData.zip -d ${TMPDIR}/Alexandre2018_SuppData
  mv ${TMPDIR}/Alexandre2018_SuppData/Ecotypes_SupplementalTables_I.xlsx ${OUTDIR}/Alexandre2018_SuppTableI.xlsx

  Rscript subscripts/extract_At_union_peaks.r ${OUTDIR}/Alexandre2018_SuppTableI.xlsx \
    | awk '$0 != ""' \
    > ${OUTDIR}/Arabidopsis_peaks_union.bed
fi


### create bed file with cis-regulatory regions of FT according to https://doi.org/10.1105/tpc.110.074682 ###
echo -e "1\t24325920\t24326300\tFT-C\n1\t24329477\t24329715\tFT-B\n1\t24331151\t24331509\tFT-A" > ${OUTDIR}/FT_regions.bed


### get photoperiod-responsive ACRs from https://doi.org/10.1093/plcell/koaa043 ###
if [ ! -f "${OUTDIR}/photoperiod_peaks.bed" ]; then
  curl "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/plcell/33/3/10.1093_plcell_koaa043/1/koaa043_supplementary_data.zip?Expires=1648657813&Signature=nBeQ1DQ2WhqHRLqXbp8FFfwhS8hDvl1dAgnRSn1Jwe~C7x-3GhO~1sOJLd4QTnzox0OLq9Hhjf3qUsfy-LdaeBXLllu9G34xC~fv850YeAdVVzNz~s8K-yhbTLns9s7clDQmB1~r8B6vv8oL0QrdSmbQ2DKDZRZYN36xm6BK~88T~utK5yzTbcK0fGC4OLWE7GlOggjtF2mstPU8ilH~ElRQepDaVRri-2tuHayVXLlwU4y1R1PCEsd7LaOrNrm5jIP5MMrmSleH4JYdMq6nQAFWZ4sIG2L6u6Ydr1BHuH28OAbhDqYPPSHOG~7b-gQ~hbQcEB13fztyYzUzOtNs5g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA" \
    > ${TMPDIR}/Tian2021_SuppData.zip
  mkdir ${TMPDIR}/Tian2021_SuppData
  unzip ${TMPDIR}/Tian2021_SuppData.zip -d ${TMPDIR}/Tian2021_SuppData
  mv ${TMPDIR}/Tian2021_SuppData/tpc.00892.2020-s02.xlsx ${OUTDIR}/Tian2021_SuppDataI.xlsx

  Rscript subscripts/extract_photoperiod_ACRs.r ${OUTDIR}/Tian2021_SuppDataI.xlsx \
    | bedops -m - \
    | awk '$0 != "" {print $0, "PP-"NR}' \
    > ${OUTDIR}/photoperiod_peaks.bed
fi


### select 170 bp regions centered on peaks ###
for FILE in "${OUTDIR}/Arabidopsis_peaks_union.bed" "${OUTDIR}/FT_regions.bed" "${OUTDIR}/photoperiod_peaks.bed"
do
  FILENAME=`basename ${FILE}`
  awk -v OFS='\t' '{mid = int(($2 + $3) / 2); print $1, mid - 85, mid + 85, $4}' ${FILE} \
    > ${TMPDIR}/${FILENAME%%_*}_ACRs.bed
done


### combine Arabidopsis union peaks, FT regions, and photoperiod ACRs (if overlapping by less than 80%) ###
sort-bed ${TMPDIR}/FT_ACRs.bed ${TMPDIR}/photoperiod_ACRs.bed \
  | bedops -n 80.1% - ${TMPDIR}/Arabidopsis_ACRs.bed \
  | sort-bed ${TMPDIR}/Arabidopsis_ACRs.bed - \
  | awk -v OFS='\t' '{if (substr($4, 1, 2) != "At") suffix = sprintf("(%s)", substr($4, 1, 2)); else suffix = ""; print $1, $2, $3, "At-"NR suffix}' \
  > ${OUTDIR}/Arabidopsis_ACRs.bed


### get ACR sequences ###
. subscripts/get_sequences.sh Arabidopsis


### get Bur-0 ACRs ###
. subscripts/get_Bur-0_ACRs.sh

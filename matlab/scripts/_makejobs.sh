#!/bin/sh
echo "Usage: makejobs job# [file#]"

PROJECT="infbench"

module purge

# Use Intel compiler
module load ${MATLAB_MODULE}
source ${HOME}/MATLAB/setpath.sh
export MATLABPATH=${MATLABPATH}
WORKDIR="${SCRATCH}/${PROJECT}"

if [ -z "${2}" ];
	then DIRID=${1};
	else DIRID=${2};
fi

FILEID=${1}
FILENAME="joblist-${FILEID}.txt"
echo "Input #: ${1}   Output file: ${FILENAME}"

VBMC22="{'beta'}"

# Default job list
PROBSET="'vbmc22'"
PROBS=${VBMC22}
DIMS="{'2D','6D','10D'}"
NOISE="'[]'"
ALGOS="{'vbmc'}"
ALGOSET="'base'"
IDS="{'1:2','3:4','5:6','7:8','9:10','11:12','13:14','15:16','17:18','19:20'}"
IDS_SINGLE="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
IDS_FIFTY="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50'}"
IDS_CENTO="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99','100'}"

##############################################################################

# Run the program
case "${1}" in
	0)      PROBS="{'lumpy'}"
		ALGOS="{'vbmc'}"
		DIMS="{'2D','4D'}"
		IDS="{'1','2'}"
		;;
	1)	PROBS="{'beta'}"
		ALGOS="{'vbmc'}"
		ALGOSET="{'basenew','probit'}"
		DIMS="{'P411','P431','P461','P491','P433'}"
		IDS=${IDS}
		;;
        2)      PROBS="{'beta'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'noisynew','noisyprobit'}"
                DIMS="{'P411','P431','P461','P491','P433'}"
		NOISE="{'me','hi'}"
                IDS=${IDS}
		;;

	11)     PROBS="{'conbananas'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'basenew','probit'}"
                DIMS="{'D2','D4','D6','D8'}"
                IDS=${IDS_SINGLE}
                ;;
	12)     PROBS="{'conbananas'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'noisynew','noisyprobit'}"
                DIMS="{'D2','D4','D6','D8'}"
		NOISE="{'me','hi'}"
                IDS=${IDS_SINGLE}
                ;;

	101)    PROBSET="{'vbmc22'}"
                PROBS="{'acerbi2012'}"
                ALGOS="{'vbmc'}"
		ALGOSET="{'basenew','probit'}"
                DIMS="{'S101'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        102)    PROBSET="{'vbmc22'}"
                PROBS="{'acerbi2012'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'noisynew','noisyprobit'}"
                DIMS="{'S1'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;

	201)    PROBSET="{'vbmc22'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'basenew','probit'}"
                DIMS="{'S101','S102'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        202)    PROBSET="{'vbmc22'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'noisynew','noisyprobit'}"
                DIMS="{'S1','S2'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;

	301)    PROBSET="{'vbmc22'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'basenew','probit'}"
                DIMS="{'S101','S102'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        302)    PROBSET="{'vbmc22'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'noisynew','noisyprobit'}"
                DIMS="{'S1','S2'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;

        401)    PROBSET="{'vbmc22'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'basenew','probit'}"
                DIMS="{'S107','S108'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        402)    PROBSET="{'vbmc22'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'noisynew','noisyprobit'}"
                DIMS="{'S107','S108'}"
                IDS=${IDS_SINGLE}
                NOISE="'me'"
                ;;

	501)    PROBSET="{'vbmc22'}"
                PROBS="{'akrami2018b'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'basenew','probit'}"
                DIMS="{'S101'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        502)    PROBSET="{'vbmc22'}"
                PROBS="{'akrami2018b'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'noisynew','noisyprobit'}"
                DIMS="{'S1'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;

	601)    PROBSET="{'vbmc22'}"
                PROBS="{'wood2010'}"
                ALGOS="{'vbmc'}"
                ALGOSET="{'noisynew','noisyprobit'}"
                DIMS="{'D1'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;



        1001)   PROBS="{'beta'}"
                ALGOS="{'inflaplace'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'P411','P431','P461','P491','P433'}"
                IDS="'1:20'"
                ;;
        1011)   PROBS="{'conbananas'}"
                ALGOS="{'inflaplace'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'D2','D4','D6','D8'}"
                IDS="'1:20'"
                ;;
        1101)   PROBSET="{'vbmc22'}"
                PROBS="{'acerbi2012'}"
                ALGOS="{'inflaplace'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101'}"
                IDS="'1:20'"
                NOISE="'[]'"
                ;;
        1201)   PROBSET="{'vbmc22'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'inflaplace'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101','S102'}"
                IDS="'1:20'"
                NOISE="'[]'"
                ;;
        1301)   PROBSET="{'vbmc22'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'inflaplace'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101','S102'}"
                IDS=${IDS}
                NOISE="'[]'"
                ;;
        1401)   PROBSET="{'vbmc22'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'inflaplace'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S107','S108'}"
                IDS=${IDS}
                NOISE="'[]'"
                ;;
        1501)   PROBSET="{'vbmc22'}"
                PROBS="{'akrami2018b'}"
                ALGOS="{'inflaplace'}"
                ALGOSET="{'base','probit'}"
                DIMS="'S101'"
                IDS=${IDS}
                NOISE="'[]'"
                ;;

        2001)   PROBS="{'beta'}"
                ALGOS="{'wsabiplus'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'P411','P431','P461','P491','P433'}"
                IDS=${IDS}
                ;;
        2011)   PROBS="{'conbananas'}"
                ALGOS="{'wsabiplus'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'D2','D4','D6','D8'}"
                IDS=${IDS_SINGLE}
                ;;
        2101)   PROBSET="{'vbmc22'}"
                PROBS="{'acerbi2012'}"
                ALGOS="{'wsabiplus'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101'}"
                IDS=${IDS}
                NOISE="'[]'"
                ;;
        2201)   PROBSET="{'vbmc22'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'wsabiplus'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101','S102'}"
                IDS=${IDS}
                NOISE="'[]'"
                ;;
        2301)   PROBSET="{'vbmc22'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'wsabiplus'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101','S102'}"
                IDS=${IDS}
                NOISE="'[]'"
                ;;
        2401)   PROBSET="{'vbmc22'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'wsabiplus'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S107','S108'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        2501)   PROBSET="{'vbmc22'}"
                PROBS="{'akrami2018b'}"
                ALGOS="{'wsabiplus'}"
                ALGOSET="{'base','probit'}"
                DIMS="'S101'"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;


        3001)   PROBS="{'beta'}"
                ALGOS="{'mfvi'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'P411','P431','P461','P491','P433'}"
                IDS=$IDS
                ;;
        3011)   PROBS="{'conbananas'}"
                ALGOS="{'mfvi'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'D2','D4','D6','D8'}"
                IDS=$IDS
                ;;
        3101)   PROBSET="{'vbmc22'}"
                PROBS="{'acerbi2012'}"
                ALGOS="{'mfvi'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        3201)   PROBSET="{'vbmc22'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'mfvi'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101','S102'}"
                IDS=${IDS}
                NOISE="'[]'"
                ;;
        3301)   PROBSET="{'vbmc22'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'mfvi'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S101','S102'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        3401)   PROBSET="{'vbmc22'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'mfvi'}"
                ALGOSET="{'base','probit'}"
                DIMS="{'S107','S108'}"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
        3501)   PROBSET="{'vbmc22'}"
                PROBS="{'akrami2018b'}"
                ALGOS="{'mfvi'}"
                ALGOSET="{'base','probit'}"
                DIMS="'S101'"
                IDS=${IDS_SINGLE}
                NOISE="'[]'"
                ;;
	
esac

##############################################################################

echo "Job items: ${PROBSET},${PROBS},${DIMS},${NOISE},${ALGOS},${ALGOSET},${IDS}"

cat<<EOF | matlab -nodisplay
addpath(genpath('${HOME}/${PROJECT}'));
currentDir=cd;
cd('${WORKDIR}');
infbench_joblist('${FILENAME}','run${DIRID}',${PROBSET},${PROBS},${DIMS},${NOISE},${ALGOS},${ALGOSET},${IDS});
cd(currentDir);
EOF

cat ${WORKDIR}/${FILENAME}
cat ${WORKDIR}/${FILENAME} | wc

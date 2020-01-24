#!/bin/sh
echo "Usage: noisyjobs job# [file#]"

PROJECT="infbench"
#source ${HOME}/MATLAB/setroot.sh

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab/2017a
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

VBMC18="{'lumpy','studentt','cigar'}"

# Default job list
PROBSET="'vbmc18'"
PROBS=${VBMC18}
DIMS="{'2D','4D','6D','8D','10D'}"
NOISE="'me'"
ALGOS="{'vbmc'}"
ALGOSET="'base'"
IDS="{'1:2','3:4','5:6','7:8','9:10','11:12','13:14','15:16','17:18','19:20'}"
IDS_SINGLE="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
IDS_FIFTY="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50'}"
IDS_CENTO="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99','100'}"

case "${1}" in
	0)      PROBS="{'lumpy'}"
		ALGOS="{'vbmc'}"
		DIMS="{'2D','4D'}"
		IDS="{'1','2'}"
		;;
	1)	ALGOSET="{'base'}"
		IDS=${IDS_SINGLE}
		;;
	2a)     ALGOSET="{'acqmi'}"
		IDS=${IDS_SINGLE}
		;;
        2b)     ALGOSET="{'acqmi2'}"
                IDS=${IDS_SINGLE}
                ;;
        2c)     ALGOSET="{'acqmi3'}"
                IDS=${IDS_SINGLE}
                ;;
        2d)     ALGOSET="{'acqmi4'}"
                IDS=${IDS_SINGLE}
                ;;
        2e)     ALGOSET="{'acqmiss'}"
                IDS=${IDS_SINGLE}
                ;;
        2f)     ALGOSET="{'acqmiss2'}"
                IDS=${IDS_SINGLE}
                ;;
        3a)     ALGOSET="{'acqopt1'}"
                IDS=${IDS_SINGLE}
                ;;
        3b)     ALGOSET="{'acqopt3'}"
                IDS=${IDS_SINGLE}
                ;;
        3c)     ALGOSET="{'acqopt10'}"
                IDS=${IDS_SINGLE}
                ;;
        3d)     ALGOSET="{'acqopt06'}"
                IDS=${IDS_SINGLE}
                ;;
      	4)      ALGOSET="{'intmeanconstacqimi'}"
		IDS=${IDS_SINGLE}
                ;;
        4b)     ALGOSET="{'acqimiqrnoise'}"
                IDS=${IDS_SINGLE}
                ;;
        4c)      ALGOSET="{'acqimiqrnoisenegquad'}"
                IDS=${IDS_SINGLE}
                ;;
        4d)     ALGOSET="{'gpthreshruncov'}"
                IDS=${IDS_SINGLE}
                ;;
        4e)     ALGOSET="{'acqf'}"
                IDS=${IDS_SINGLE}
                ;;
        4f)     ALGOSET="{'acqfv2'}"
                IDS=${IDS_SINGLE}
                ;;
        4g)     ALGOSET="{'acqf2reg'}"
                IDS=${IDS_SINGLE}
                ;;
        4h)     ALGOSET="{'acqfreg2','acqpropreg2'}"
                DIMS="{'20D'}"
                IDS=${IDS_SINGLE}
                ;;
        4h1)    ALGOSET="{'K1'}"
                IDS=${IDS_SINGLE}
		;;
	4h2)    ALGOSET="{'K2'}"
                IDS=${IDS_SINGLE}
                ;;
        4h5)    ALGOSET="{'K5'}"
                IDS=${IDS_SINGLE}
                ;;
        4i)     ALGOSET="{'oldsettings'}"
                IDS=${IDS_SINGLE}
                ;;
        4j)     ALGOSET="{'step1'}"
                IDS=${IDS_SINGLE}
                ;;
        4k)     ALGOSET="{'step5'}"
                IDS=${IDS_SINGLE}
                ;;
        4l)     ALGOSET="{'acqmistep5'}"
                IDS=${IDS_SINGLE}
                ;;
        4m)     ALGOSET="{'acq2step1'}"
                IDS=${IDS_SINGLE}
                ;;
        4n)     ALGOSET="{'se'}"
                IDS=${IDS_SINGLE}
                ;;
        4o)     ALGOSET="{'test'}"
                IDS=${IDS_SINGLE}
                ;;
       	4p)    	ALGOSET="{'finalK'}"
                IDS=${IDS_SINGLE}
                ;;
        4q)     ALGOSET="{'step1'}"
                IDS=${IDS_SINGLE}
                ;;
        4r)     ALGOSET="{'newdef2'}"
                IDS=${IDS_SINGLE}
                ;;
        4s)     ALGOSET="{'step1mi'}"
                IDS=${IDS_SINGLE}
                ;;
        4t)     ALGOSET="{'newdef3'}"
                IDS=${IDS_SINGLE}
                ;;
        4u)     ALGOSET="{'outnoise'}"
                IDS=${IDS_SINGLE}
                ;;
        5)      ALGOSET="{'noempgp'}"
                IDS=${IDS_SINGLE}
		;;
        6)      ALGOSET="{'noempgp'}"
                IDS=${IDS_SINGLE}
		;;
	7)	ALGOSET="{'morelbomid'}"
		IDS=${IDS_SINGLE}
		;;
        8)      ALGOSET="{'acqprop'}"
                ;;
        9)      ALGOSET="{'acqpropnorot'}"
		IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
        10)     ALGOSET="{'acqpropfewrot'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
        11)     ALGOSET="{'acqfnorot'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
 		;;
        12)     ALGOSET="{'acqpropfnorot'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
       13)     ALGOSET="{'acqpropf'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
       14)     ALGOSET="{'acqpropfnorotbz'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
       15)     ALGOSET="{'acqpropfbz'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
       16)     ALGOSET="{'adaptive'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
       17)     ALGOSET="{'adaptiveless'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
       18)     ALGOSET="{'adaptivelessnogp'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
       19)     ALGOSET="{'basenogp'}"
                IDS="{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}"
                ;;
       20)     ALGOSET="{'basekone'}"
               IDS=$IDS_SINGLE
                ;;
       21)     	ALGOSET="{'acqproponly'}"
               	IDS=$IDS_SINGLE
		;;
       22)      ALGOSET="{'acqus'}"
                IDS=$IDS_SINGLE
                ;;
       23)      ALGOSET="{'acqev'}"
                IDS=$IDS_SINGLE
                ;;
       24)      ALGOSET="{'detvars3'}"
                IDS=$IDS_SINGLE
                ;;
       25)      ALGOSET="{'detvars5'}"
                IDS=$IDS_SINGLE
                ;;
       26)      ALGOSET="{'acqpropcontrol'}"
                IDS=$IDS_SINGLE
                ;;

	30)     ALGOS="{'parallelgp'}"
                IDS=${IDS_SINGLE}
                ;;
        30b)    ALGOS="{'parallelgp'}"
                IDS=${IDS_SINGLE}
		NOISE="'hi'"
                ;;

        31)     ALGOS="{'parallelgp'}"
		ALGOSET="{'fast'}"
                IDS=${IDS_SINGLE}
                ;;


	50)     ALGOS="{'wsabi'}"
		IDS=${IDS_SINGLE} 
                ;;
        51)     ALGOS="{'wsabi'}"
		ALGOSET="{'mm'}"
		IDS=${IDS_SINGLE}
                ;;
        52)     ALGOS="{'wsabi'}"
		ALGOSET="{'search'}"
                IDS=${IDS_SINGLE}
                ;;
        53)     ALGOS="{'wsabi'}"
                ALGOSET="{'mm2'}"
                IDS=${IDS_SINGLE}
                ;;
        55)     ALGOS="{'bbq'}"
		IDS=${IDS_SINGLE}
                ;;
        56)     ALGOS="{'bbq'}"
                ALGOSET="{'marginal'}"
		IDS=${IDS_SINGLE}
                ;;
	60)     ALGOS="{'bmc'}"
                ;;
        70)     ALGOS="{'smc'}"
                ;;
        80)     ALGOS="{'ais'}"
                ;;
	90)	ALGOS="{'agp'}"
		IDS=$IDS_SINGLE
		;;
       	91)     ALGOS="{'agp@long'}"
                IDS=$IDS_SINGLE
                ;;
        92)     ALGOS="{'agp@reg2'}"
                IDS=$IDS_SINGLE
                ;;
        95)     ALGOS="{'bape'}"
                IDS=$IDS_SINGLE
                ;;
        96)     ALGOS="{'bape@step1'}"
                IDS=$IDS_SINGLE
                ;;
        97)     ALGOS="{'bape@nqreg2'}"
                IDS=$IDS_SINGLE
                ;;
        98)     ALGOS="{'bape@reg'}"
                IDS=$IDS_SINGLE
                ;;
        101)    PROBS="{'goris2015'}"
		ALGOS="{'laplace'}"
                DIMS="{'S7','S8','S9','S10','S11','S12'}"
                IDS="{'1','2','3'}"
                ;;
        102)    PROBS="{'goris2015'}"
                ALGOS="{'wsabi','wsabi@mm','bmc','smc','ais'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
	102b)   PROBS="{'goris2015'}"
                ALGOS="{'bbq','bbq@marginal'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        102c)   PROBS="{'goris2015'}"
                ALGOS="{'wsabi','wsabi@mm'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        103)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@base'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        103L)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@baseluigi'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        104)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqus'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104b)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqev'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104c)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqi'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104d)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqf1'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104e)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqf'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104f)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqfv2'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104g)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqf2reg'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104h)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@K1','vbmc@K2','vbmc@K5'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        104i)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@oldsettings'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104j)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@step1'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104k)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@step5'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104l)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqmistep5'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104m)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acq2step1'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104n)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@se'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104o)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@test'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104p)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@finalK'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104q)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@step1'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104r)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@newdef2'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104s)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@step1mi'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104t)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@newdef3'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104u)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@outnoise'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
	104v)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@band'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        104w)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@fixgpmean'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        105)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@control'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        106)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@noempgp'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        107)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@morelbomid'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        108)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@cheapgpmid'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        109)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqproponly'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        110)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqusonly'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        111)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqvusonly'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        112)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqpropnoprune'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        113)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqusnoprune'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        114)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@acqpropcontrol'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        115)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@cheapgpmidlessprune'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
	116)    PROBS="{'goris2015'}"
                ALGOS="{'vbmc@gpthreshmidlessprune'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        116b)   PROBS="{'goris2015'}"
                ALGOS="{'vbmc@lessprune'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        120)    PROBS="{'goris2015'}"
                ALGOS="{'agp'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        121)    PROBS="{'goris2015'}"
                ALGOS="{'agp@long'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        122)    PROBS="{'goris2015'}"
                ALGOS="{'agp@reg2'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
	130)    PROBS="{'goris2015'}"
                ALGOS="{'bape'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        131)    PROBS="{'goris2015'}"
                ALGOS="{'bape@step1'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_CENTO
                ;;
        132)    PROBS="{'goris2015'}"
                ALGOS="{'bape@nqreg2'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        133)    PROBS="{'goris2015'}"
                ALGOS="{'bape@reg'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;
        140)    PROBS="{'goris2015'}"
                ALGOS="{'wsabi@search'}"
                DIMS="{'S7','S8'}"
                IDS=$IDS_FIFTY
                ;;


	201)    PROBS="{'yacht'}"
                ALGOS="{'laplace'}"
                DIMS="{'8D'}"
                IDS="{'1','2','3'}"
                ;;
        202)    PROBS="{'yacht'}"
                ALGOS="{'wsabi','wsabi@mm','bmc','smc','ais'}"
                DIMS="{'8D'}"
                IDS=$IDS_FIFTY
                ;;
        203)    PROBS="{'yacht'}"
                ALGOS="{'bbq','bbq@marginal'}" #,'agp','bape','bape@negquad'}"
                DIMS="{'8D'}"
                IDS=$IDS_FIFTY
                ;;
        204)    PROBS="{'yacht'}"
                ALGOS="{'vbmc','vbmc@acqfreg2','vbmc@acqpropreg2'}"
                DIMS="{'8D'}"
                IDS=$IDS_CENTO
                ;;
        204b)   PROBS="{'yacht'}"
                ALGOS="{'vbmc@oldsettings','vbmc@step1','vbmc@acqmistep1','vbmc@step5','vbmc@acqmistep5'}"
                DIMS="{'8D'}"
                IDS=$IDS_CENTO
                ;;
	205)    PROBS="{'yacht'}"
                ALGOS="{'wsabi@search'}"
                DIMS="{'8D'}"
                IDS=$IDS_CENTO
                ;;
        206)    PROBS="{'yacht'}"
                ALGOS="{'bape@step1'}"
                DIMS="{'8D'}"
                IDS=$IDS_CENTO
                ;;



        304b)   PROBS="{'funnel'}"
                ALGOS="{'vbmc@oldsettings','vbmc@step1','vbmc@acqmistep1','vbmc@step5','vbmc@acqmistep5'}"
                IDS=${IDS_SINGLE}
		;;
        305)    PROBS="{'funnel'}"
                ALGOS="{'wsabi@search'}"
                IDS=${IDS_SINGLE}
		;;
        306)    PROBS="{'funnel'}"
                ALGOS="{'bape@step1'}"
                #IDS=${IDS_SINGLE}
                IDS="{'1','2'}"
		;;



esac

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

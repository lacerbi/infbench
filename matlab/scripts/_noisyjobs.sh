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
DIMS="{'2D','6D','10D'}"
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
        4b)     ALGOSET="{'rotoup'}"
                IDS=${IDS_SINGLE}
                ;;
        4b2)    ALGOSET="{'rotoup'}"
                IDS=${IDS_SINGLE}
		NOISE="{'hi'}"
                ;;
        4c)     ALGOSET="{'acqimiqrnoiserotonomcmc2'}"
                IDS=${IDS_SINGLE}
                ;;
        4c2)    ALGOSET="{'acqimiqrnoiserotonomcmc2'}"
                IDS=${IDS_SINGLE}
		NOISE="{'hi'}"
                ;;
        4d)     ALGOSET="{'acqnoise'}"
                IDS=${IDS_SINGLE}
                ;;
        4d2)    ALGOSET="{'acqnoise'}"
                IDS=${IDS_SINGLE}
                NOISE="{'hi'}"
                ;;
        4e)     ALGOSET="{'acqimiqrnoiseup'}"
                IDS=${IDS_SINGLE}
                ;;
        4e2)    ALGOSET="{'acqimiqrnoiseup'}"
                IDS=${IDS_SINGLE}
                NOISE="{'hi'}"
                ;;
        4f)     ALGOSET="{'acqimiqrnoiserotoupthin100'}"
                IDS=${IDS_SINGLE}
                ;;
        4f2)    ALGOSET="{'acqimiqrnoiserotoupthin100'}"
                IDS=${IDS_SINGLE}
                NOISE="{'hi'}"
                ;;
        5a)     ALGOSET="{'oldsettings'}"
                IDS=${IDS_SINGLE}
		;;
        5a2)    ALGOSET="{'oldsettings'}"
                IDS=${IDS_SINGLE}
                NOISE="{'hi'}"
                ;;
        6a)     ALGOSET="{'newbase'}"
                IDS=${IDS_SINGLE}
		;;
        6a2)     ALGOSET="{'newbase'}"
                IDS=${IDS_SINGLE}
		NOISE="{'hi'}"
                ;;
	7a)	ALGOSET="{'renewbase2d'}"
		IDS=${IDS_SINGLE}
		;;
        7a2)     ALGOSET="{'renewbase2d'}"
                IDS=${IDS_SINGLE}
		NOISE="{'hi'}"
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

        404)    PROBSET="{'vbmc20'}"
		PROBS="{'wood2010'}"
                ALGOS="{'vbmc@newbase'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
		NOISE="'[]'"
                ;;
        405)    PROBSET="{'vbmc20'}"
                PROBS="{'wood2010'}"
                ALGOS="{'vbmc@renewbase2d'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        405b)    PROBSET="{'vbmc20'}"
                PROBS="{'wood2010'}"
                ALGOS="{'vbmc@acqimiqrnoisewarp'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        406)    PROBSET="{'vbmc20'}"
                PROBS="{'wood2010'}"
                ALGOS="{'vbmc@rotoup'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        407)    PROBSET="{'vbmc20'}"
                PROBS="{'wood2010'}"
                ALGOS="{'vbmc@acqimiqrnoiserotocorrupthin100'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        408)    PROBSET="{'vbmc20'}"
                PROBS="{'wood2010'}"
                ALGOS="{'vbmc@oldsettings'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        430)    PROBSET="{'vbmc20'}"
                PROBS="{'wood2010'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        504)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@newbase'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        505)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@renewbase2d'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        505b)   PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@acqimiqrnoisewarp'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        506)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@rotoup'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        507)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@acqimiqrnoiserotocorrupthin100'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        508)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@oldsettings'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        530)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
	550)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@newdef'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        551)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@roto'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
	552)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@acqboth'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        553)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'vbmc@acqbothup'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
	560)    PROBSET="{'vbmc20'}"
                PROBS="{'acerbidokka2018'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;

        
        604)    PROBSET="{'vbmc20'}"
                PROBS="{'price2018'}"
                ALGOS="{'vbmc@newbase'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        605)    PROBSET="{'vbmc20'}"
                PROBS="{'price2018'}"
                ALGOS="{'vbmc@renewbase2d'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        606)    PROBSET="{'vbmc20'}"
                PROBS="{'price2018'}"
                ALGOS="{'vbmc@rotoup'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        607)    PROBSET="{'vbmc20'}"
                PROBS="{'price2018'}"
                ALGOS="{'vbmc@acqimiqrnoiserotocorrupthin100'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        608)    PROBSET="{'vbmc20'}"
                PROBS="{'price2018'}"
                ALGOS="{'vbmc@oldsettings'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        608b)   PROBSET="{'vbmc20'}"
                PROBS="{'price2018'}"
                ALGOS="{'vbmc@acqimiqrnoiserotoupthin100'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        630)    PROBSET="{'vbmc20'}"
                PROBS="{'price2018'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;

        704)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@newbase'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        705)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@renewbase2d'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        707)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@acqimiqrnoiserotocorrupthin100'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        708)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@oldsettings'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        708b)   PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@acqimiqrnoiserotoupthin100'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        730)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'1','2'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        750)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@newdef'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        751)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@roto'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        752)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@acqboth'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        753)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'vbmc@acqbothup'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        760)    PROBSET="{'vbmc20'}"
                PROBS="{'krajbich2010'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'101','102'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;

        804)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@newbase'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        805)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@renewbase2d'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        807)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@acqimiqrnoiserotocorrupthin100'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        807b)   PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@acqimiqrnoiseup2'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        808)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@oldsettings'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        808b)   PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@acqnoisemcmc'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        808c)   PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@acqimiqrnoiserotoup3'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        808d)   PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@acqimiqrnoiserotoupthin100'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        808e)   PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@acqimiqrnoiserotoup5'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        830)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'7','8'}"
                IDS=$IDS_CENTO
                NOISE="'me'"
                ;;
        850)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@newdef'}"
                DIMS="{'107','108'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        851)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@roto'}"
                DIMS="{'107','108'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        852)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@acqboth'}"
                DIMS="{'107','108'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        853)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'vbmc@acqbothup'}"
                DIMS="{'107','108'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        860)    PROBSET="{'vbmc20'}"
                PROBS="{'goris2015b'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'107','108'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;


        904)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@newbase'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        905)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@renewbase2d'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        907)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@acqimiqrnoiserotocorrupthin100'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        908)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@oldsettings'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        908d)   PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@acqimiqrnoiserotoupthin100'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        908e)   PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@acqimiqrnoiserotoup5'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        930)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'1'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        950)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@newdef'}"
                DIMS="{'101'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        951)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@roto'}"
                DIMS="{'101'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        952)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@acqboth'}"
                DIMS="{'101'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        953)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'vbmc@acqbothup'}"
                DIMS="{'101'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
                ;;
        960)    PROBSET="{'vbmc20'}"
                PROBS="{'akrami2018'}"
                ALGOS="{'parallelgp'}"
                DIMS="{'101'}"
                IDS=$IDS_CENTO
                NOISE="'[]'"
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

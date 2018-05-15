{\rtf1\ansi\ansicpg1252\cocoartf1561\cocoasubrtf200
{\fonttbl\f0\fmodern\fcharset0 Courier;\f1\fmodern\fcharset0 Courier-Bold;\f2\fnil\fcharset0 Menlo-Regular;
\f3\froman\fcharset0 TimesNewRomanPSMT;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red0\green0\blue0;\red255\green255\blue255;
}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\csgray\c0;\csgray\c100000;
}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\sl280\partightenfactor0

\f0\fs24 \cf2 \expnd0\expndtw0\kerning0
% 
\f1\b make
\f0\b0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardeftab720\pardirnatural\partightenfactor0

\f2\fs22 \cf3 \cb4 \kerning1\expnd0\expndtw0 \CocoaLigature0 icpc -I$\{TACC_PAPI_INC\} -O0 scheduler.cpp $\{TACC_PAPI_LIB\}/libpapi.a -o para1 -std=c++11 -ltbb\
icpc -I$\{TACC_PAPI_INC\} -O0 scheduler.cpp $\{TACC_PAPI_LIB\}/libpapi.a -o para3 -std=c++11 -ltbb\
sbatch task1.script\
sbatch task3.script\
\
\pard\pardeftab720\partightenfactor0

\f3\i \cf0 \cb1 \expnd0\expndtw0\kerning0
\CocoaLigature1 icpc -I/opt/papi/intel/include -O0 part2.cpp /opt/papi/intel/lib/libpapi.a -o part1.o -ltbb -lrt}
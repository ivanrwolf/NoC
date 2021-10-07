#! /bin/bash

# Check if input file is provided
if [ $# -eq 0 ]
    then
    echo "No arguments supplied"
    echo "Please provide and network file in ncol format"
    exit
fi

# Check weka default location and run
if [ -f "./weka-3-8-5/weka.jar" ]; then
    echo "./weka-3-8-5/weka.jar exists"
else 
    echo "./weka-3-8-5/weka.jar does not exist"
    exit
fi

# Paths
path_input=$1
path_weka_model="Final_model.model"
path_weka_input="04_file_for_classification.arff"
path_test_predictions="05_predictions.tsv"

# Run discretizer
echo "Calculating topologies and discretizing"
wait
python grn.py ${path_input}
wait

# Run weka
echo "Making predictions"
wait
java -Xmx15200M -cp ./weka-3-8-5/weka.jar weka.classifiers.meta.Vote -l ${path_weka_model} -p 0 -T ${path_weka_input} | grep -v '^$' |sed '1d' | sed -e 's/1\://g' -e 's/2\://g' -e 's/[space]*/\t/'| awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' > ${path_test_predictions}.temp
wait

echo "Reformating output"
wait
echo -e "gene	inst	actual	predicted	error	prediction" > ${path_test_predictions}
wait
cat ${path_test_predictions}.temp | grep -w "+" >> ${path_test_predictions}.tmp2
wait
cat ${path_test_predictions}.temp | grep -vw "+" | grep -v "error" | awk '{OFS="\t"}{print $1,$2,$3, "-",$4}' >> ${path_test_predictions}.tmp2
wait
cat 03_network_parameters_discretized.csv | cut -d "," -f1 | sed "1d" > ${path_test_predictions}.genes.tmp
wait
paste ${path_test_predictions}.genes.tmp ${path_test_predictions}.tmp2 | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' >> ${path_test_predictions}
wait


echo "Cleaning up temporary files"
wait
rm ${path_test_predictions}.temp ${path_test_predictions}.tmp2 ${path_test_predictions}.genes.tmp
wait

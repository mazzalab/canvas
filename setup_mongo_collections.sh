#!/bin/bash
usage() { echo "Usage: $0 [-d <string>] [-r <string>]" 1>&2; exit 1; }

die () {
    echo "ERROR: $*. Aborting." >&2
    exit 1
}

dir=''
build=false
cleanup=false
regex1="genelist"
regex2="tarbase"
regex3="targetscan"
regex4="mirbase"
regex5="enhancer"
regex6="overlap"
while getopts "d:r:" opt; do
      case $opt in
            d) $cleanup && die "Cannot specify option -d with -r"
               build=true
               dir=${OPTARG};;

            r) $build && die "Cannot specify option -d with -r"
               cleanup=true
               db=${OPTARG};;
            \?)
                usage
                ;;
      esac
    shift $((OPTIND-1))

    if [[ "${cleanup}" == true ]] ; then
        mongo ${db} --eval "db.dropDatabase()"


    elif [[ "${build}" == true ]] ; then
        echo "build"
        echo "importing from $dir"
        echo "db will be named $(basename $dir)"
        for file in ${dir}/*.txt; do
            s=${file##*/};
            echo "Importing $file as ${s%.txt}";
            if [[ "$file" =~ $regex1 ]]; then
                mongoimport --host=127.0.0.1 -d $(basename $dir) -c ${s%_genelist.txt} --type tsv --file $file --fields gene
            elif [[ "$file" =~ $regex2 ]]; then
                mongoimport --host=127.0.0.1 -d $(basename $dir) -c ${s%.txt} --type tsv --file $file --fields geneID,geneName,mirna
            elif [[ "$file" =~ $regex3 ]]; then
                mongoimport --host=127.0.0.1 -d $(basename $dir) -c ${s%.txt} --type tsv --file $file --fields geneName,mirna
            elif [[ "$file" =~ $regex4 ]]; then
                mongoimport --host=127.0.0.1 -d $(basename $dir) -c ${s%.txt} --type tsv --file $file --fields chr,source,feature,start,end,score,strand,frame,info
            elif [[ "$file" =~ $regex5 ]]; then
                mongoimport --host=127.0.0.1 -d $(basename $dir) -c ${s%.txt} --type tsv --file $file --fields info,chr,start,end,CellType,source,TargetGene,ConservationMean,ConservationMedian
                mongo $(basename $dir) --eval "db.${s%.txt}.createIndex({chr:1});"
            elif [[ "$file" =~ $regex6 ]]; then
                mongoimport --host=127.0.0.1 -d $(basename $dir) -c ${s%.txt} --type tsv --file $file --fields chr,start,end,type
                mongo $(basename $dir) --eval "db.${s%.txt}.createIndex({chr:1});"
            else
                mongoimport --host=127.0.0.1 -d $(basename $dir) -c ${s%.txt} --type tsv --file $file --fields chr,start,end,info
                mongo $(basename $dir) --eval "db.${s%.txt}.createIndex({chr:1});"
            fi

        done
    fi

done

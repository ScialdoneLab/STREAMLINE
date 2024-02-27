#!/bin/bash

experimental='false'
directed=''
local='false'
global='false'
statistical='false'


while getopts 'edslg' flag; do
    case "${flag}" in
        e) experimental='true' ;;
        d) directed='--directed' ;;
        s) statistical='true' ;;
        l) local='true' ;;
        g) global='true' ;;
    esac
done

echo "Running STREAMLINE..."
cd STREAMLINE/

if [ $global ==  "true" ]; then
    echo "Evaluating synthetic networks (global)"
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_curated.yaml --STREAMLINE_global $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_ER.yaml --STREAMLINE_global $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SF.yaml --STREAMLINE_global $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SSF.yaml --STREAMLINE_global $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SW.yaml --STREAMLINE_global $directed
fi

if [ $local ==  "true" ]; then
    echo "Evaluating synthetic networks (local)"
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_curated.yaml --STREAMLINE_local $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_ER.yaml --STREAMLINE_local $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SF.yaml --STREAMLINE_local $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SSF.yaml --STREAMLINE_local $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SW.yaml --STREAMLINE_local $directed
fi

if [ $statistical ==  "true" ]; then
    echo "Evaluating synthetic networks (EPr)"
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_curated.yaml -e $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_ER.yaml -e $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SF.yaml -e $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SSF.yaml -e $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SW.yaml -e $directed
    
    echo "Evaluating synthetic networks (AUC)"
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_curated.yaml -a $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_ER.yaml -a $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SF.yaml -a $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SSF.yaml -a $directed
    python  BLEvaluator_adapted.py --config config-files/synthetic/config_SW.yaml -a $directed
fi

if [ $experimental ==  "true" ]; then
    if [ $global ==  "true" ]; then
        echo "Evaluating experimental networks (global)"
        python  BLEvaluator_adapted.py --config config-files/experimental/config_yeast.yaml --STREAMLINE_global $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_hESC.yaml --STREAMLINE_global $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_mDC.yaml --STREAMLINE_global $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_mESC.yaml --STREAMLINE_global $directed
    fi
    
    if [ $local ==  "true" ]; then
        echo "Evaluating experimental networks (local)"
        python  BLEvaluator_adapted.py --config config-files/experimental/config_yeast.yaml --STREAMLINE_local $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_hESC.yaml --STREAMLINE_local $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_mDC.yaml --STREAMLINE_local $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_mESC.yaml --STREAMLINE_local $directed
    fi
    
    if [ $statistical ==  "true" ]; then
        echo "Evaluating experimental networks (EPr)"
        python  BLEvaluator_adapted.py --config config-files/experimental/config_yeast.yaml -e $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_hESC.yaml -e $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_mDC.yaml -e $directed
        python  BLEvaluator_adapted.py --config config-files/experimental/config_mESC.yaml -e $directed
    fi
fi





#!/bin/bash
if test -e log.lammps 
	then 
                echo "Start scraping for relevent data."

		#########################################################
		# scrape the data
		grep "Temp" log.lammps | tr "T" "\n" | grep emp | replace emp Temp | replace "=" "" > log.csv
		
		# test if the data is indeed scraped
		if test -n $(grep "Temp" log.csv | wc -l )
		then
		  echo "Scraped temperature data."
		else 	
		  echo "Temperature data not scraped."
		fi 
		##########################################################
		grep Volume log.lammps | tr "=" " " >> log.csv	
		
		# test if the data is indeed scraped
		if test -n $(grep "Volume" log.csv | wc -l )
		then
		  echo "Scraped volume data."
		else 	
		  echo "Volume data not scraped."
		fi 
		##########################################################
		grep "Press"  log.lammps | tr "P" "\n" | grep ress | replace ress Press | replace "=" "" >> log.csv
		# test if the data is indeed scraped
		if test -n $(grep "Press" log.csv | wc -l )
		then
		  echo "Scraped pressure data."
		else 	
		  echo "Pressure data not scraped."
		fi 		
		##########################################################
		grep "TotEng" log.lammps | tr "K" "\n" | grep "TotEng" | replace "=" "" >> log.csv
		# test if the data is indeed scraped
		if test -n $(grep "TotEng" log.csv | wc -l )
		then
		  echo "Scraped total energy data."
		else 	
		  echo "Total energy data not scraped."
		fi 
		##########################################################
		grep "PotEng" log.lammps | replace "E_bond" "This" | tr "This" "\n" | grep "PotEng" | replace "=" "" >>log.csv 
		
		# test if the data is indeed scraped
		if test -n $(grep "PotEng" log.csv | wc -l )
		then
		  echo "Scraped potential data."
		else 	
		  echo "Potential data not scraped."
		fi 
		##########################################################
		grep KinEng log.lammps | tr "Temp" "\n" | tr "K" "\n" | grep inEng | replace "inEng" "KinEng" | replace "=" "" >> log.csv 
		# test if the data is indeed scraped
		if test -n $(grep "KinEng" log.csv | wc -l )
		then
		  echo "Scraped kinetic energy data."
		else 	
		  echo "Kinetic energy data not scraped."
		fi 
		##########################################################
		echo "The total time used is $SECONDS."
		echo "The size of the scraped file is $(du -sh log.csv | tr "log.csv" "\n" | grep K)" 
	else 
		echo "The file log.lammps is not here."
fi

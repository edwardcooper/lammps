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
		echo "Finish scraping volume."

		grep "Press"  log.lammps | tr "P" "\n" | grep ress | replace ress Press | replace "=" "" >> log.csv
                echo "Finish scraping pressure." 
		
		grep "TotEng" log.lammps | tr "K" "\n" | grep "TotEng" | replace "=" "" >> log.csv
		echo "Finish scraping total energy"
		
		grep "PotEng" log.lammps | replace "E_bond" "This" | tr "This" "\n" | grep "PotEng" | replace "=" "" >>log.csv 
		echo "Finish scraping potential energy."

		grep KinEng log.lammps | tr "Temp" "\n" | tr "K" "\n" | grep inEng | replace "inEng" "KinEng" | replace "=" "" >> log.csv 
		echo "Finish scraping kinetic energy."
	else 
		echo "The file log.lammps is not here."
fi


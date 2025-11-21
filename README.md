# [HLA-SeqPredict](https://ext5.nbt.or.th)

> A web application designed to streamline the process for ADR risk assesment

## About 

This repository contains 3 directories that are a part of my HLA-SeqPredict project.

* `ADRDatabaseSetup` Handles any python files that were used to create / handle any internal databases containing information on ADR.
* `WebFrontend` The code for the frontend of the website. Written using the Vue.js framework.
* `WebBackend` The code for the backend of the website. Written using the Flask framework.

This project is part of the [Satit Kaset Senior Project Program](https://academic.kus.ku.ac.th/seniorproject/searchByGrade.php?command=SEARCH&txtKeyword=HLA-SeqPredict).

## Frontend Website

The frontend website can be viewed at https://ext5.nbt.or.th. This repository will contain the source code of that website.

### How it works

This frontend website acts as a GUI for the user to directly interact with. It directly communicates  with the backend of the website using http requests that send user provided information and the necessary options to run tools hosted in the backend. It also serves as a nice interface to visualize data as well.

### Web Pages

The website currently has 8 pages.
 1. **Home page** This contains the home page of the website and some statistics trivia.
![Home page](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/home.png)

 2. **Login page** This page contains a page for the user to login.
 ![Login page](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/login.png)
 
 3. **Signup page** This page allows the user to signup with a new account. The user must provide an identifiable email (There is no method for verification, it only serves as a unique identifier like a username.) and a display name.
 ![Login page](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/signup.png)

4. **Profile page** This page allows the user to modify account information.
![Profile page](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/profile.png)
 
 5. **Upload page** This page handles file uploads. It can handle FastQ DNA sequence files and can handle both single end and paired end data. All files are securely stored on the server until the user wishes to delete it.
 ![Upload  page](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/upload.png)
 
 6. **Dashboard page** This page will show the user all the current stored sequences and what their status is. The user can perform many actions on the sequences they uploaded. These actions include viewing the logs of the commands that are ran in the backend, visualizing with IGV, running a [Haplotyper](https://github.com/DiltheyLab/HLA-LA), getting the report in the form of a pdf, and permanently deleting the file of the servers.
 ![Dashboard page](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/dashboard.png)
 
 7. **IGV page** This page is the redirected page once the user presses `Visuzlize with IGV` on the Dashboard page. This page contains the web embedded IGV.js program that will automatically visualize the uploaded file to the 6th human chromosome (chr6).
 ![IGV page](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/igv.png)
 
 8. **Search page**  This page contains 2 graphs. A histogram that shows the frequencies of people having a set haplotype. It can be filtered based on the region (Which will show what that regions percent chance is) and the related drug correlation (If it relates to having an ADR to that drug). 
 ![Histogram](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/histogram.png)
	 The other graph is a sunburst chart that displays information in a level-like manner. These "levels" can filter between Disease, Drug, ADR Gene Loci, and Allele.
 ![Sunburst](https://raw.githubusercontent.com/theme222/HLA-SeqPredict/main/Pictures/sunburst.png)
 
## Backend website 

This backend website houses all the internal tools and logic required for the services provided by the frontend. 

### How it works

The backend is made using python with a Flask backend framework to handle communication. There are also 3 external programs that is installed in the backend servers: [BWA](https://github.com/lh3/bwa) [Samtools](https://www.htslib.org/) and [HLA*LA](https://github.com/DiltheyLab/HLA-LA).

### API End points
* "/" Root home 
  * check if the server is available
* "/api/signup" `POST`
  * *required data* : name email password
  * Signs up user to database
* "/api/login" `POST` 
  * *required data* : email password
  * Logs user into website and gives a cookie :) 
  * example: `de8f06ba571969ae25baff4f39dbaf70__1734145084`
* "/api/checkEmail" `POST` 
  * *required data* : email
  * Checks if email is valid (email must be unique)
* "/api/getAccountInfo" `POST`
  * *required data* : cookie
  * Gives account info (registered name, email, id)
* "/api/chanceAccountInfo" `POST`
  * *required data* : cookie, *[name, email, password] (optional)*
  * Changes account info based on data given
* "/api/getSequences" `POST`
  * *required data* : cookie
  * Lists sequence information
* "/api/uploadFile" `POST`
  *  *required data* : file, sequence label,cookie 
  * Uploads file to database
* "/api/run/hla_la" `POST`
  * *required data* : cookie, sequence label
  * Runs HLA*LA on sequence
* "/api/getResults/hla_la" `POST`
  * *required data* : cookie, sequence label
  * Gets typing results from HLA*LA (If completed)
* "/api/requestFile/igv" `POST`
  * *required data* : sequence label, cookie
  * Opens up a channel for loading files required by IGV.js
  * Files : xxx.bam, xxx.bam.bai
* "/api/deleteSequence" `POST`
  * *required data* : cookie, sequence label
  * Deletes selected sequence
* "/data/sunburst" `POST`
  * *required data* : country_filter, minP_filter
  * Returns sunburst json object based on provided filters. 
  * ***LONG RUNTIME (20 ms) AND LARGE DATA (2 MB)***
* "/statistic" `GET`
  *  Retrieves stats from the backend
  * Total number of files uploaded
* "/download/~token~" `GET`
  * *~token~* : the token you get from "/api/requestFile/igv"
* "/reference/chr6.fa" `GET`
  * Downloads reference 
* "/reference/chr6.fa.fai" `GET`
  * Downloads reference


# Back end website made with Flask
## Info
This is a backend website made using the Flask package
it is not intended to be run individually and should be paired with the frontend.
---
## API Channels
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
  * [Total number of files uploaded]
* "/download/<token>" `GET`
  * *<token>* : the token you get from "/api/requestFile/igv"
* "/reference/chr6.fa" `GET`
  * Downloads reference 
* "/reference/chr6.fa.fai" `GET`
  * Downloads reference

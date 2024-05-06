<script setup>
/* eslint-disable */
import { ref } from 'vue'
import { reactive } from 'vue'
import axios from 'axios'
import Cookies from 'js-cookie'
import { useRouter, useRoute } from 'vue-router'
import { getAccountInfo } from '@/scripts/AccountFunc'

const router = useRouter()
const route = useRoute()

const fileName = ref(null)
const fileSize = ref("0 Kb")
const accountName = ref(null)
const validLabel = ref(true)
const labelValue = ref('')

function changeFileInfo(event){
  fileName.value = event.target.files[0].name
  fileSize.value = formatFileSize(event.target.files[0].size)
}

const formatFileSize = (bytes) => {
  const units = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
  let index = 0;
  while (bytes >= 1024 && index < 4) {
    bytes /= 1024;
    index++;
  }
  return `${bytes.toFixed(0)} ${units[index]}`
};

function checkValidLabel(){
  // Define the list of characters to check for
  const charList = ['!', '@', '#', '$', '%', '^', '&', '*', '(', ')'];

  // Create a regular expression pattern to match any character in the list
  const pattern = new RegExp('[' + charList.join('') + ']');

}

function uploadFile(){
  console.log(document.getElementById('file').files[0])
    if (document.getElementById('file').files[0] === undefined){
      alert("No file specified") 
      return 0
    }
    const formData = new FormData()
    let filenameToSet  = String(Date.now()/1000)
    if (labelValue.value !== "") {filenameToSet = labelValue.value}
    formData.append('file', document.getElementById('file').files[0], filenameToSet)
    axios.post("http://localhost:7000/api/uploadFile", 
              formData, 
              {headers: {"Content-Type" : 'multipart/form-data', 
                        'Authorization': `Bearer ${Cookies.get('session_cookie')}`} })
    .then(response => {
      alert(response.data)
      router.push({name:'dashboard'})
    })
    .catch(error => {
      alert("Please use alpha numeric charectors and _ only")
      console.error('Error uploading file  :', error);
    });
  }


getAccountInfo().then(name => {accountName.value = name})
</script>

<template>

  <div class="background"></div>

  <div class="titleWrapper">
      
      <h1 class="titleText">
          <RouterLink to="home">HLA-ADR Prediction</RouterLink>
      </h1>
      
      <label class="accountName">{{accountName || "No account"}}</label>

      
  </div>
  
  <div class="headingWrapper">
      <h2 class="heading">
        Upload your DNA file (.FastQ)
      </h2>
  </div>
  
  <div class="uploadButtonWrapper">
    <div style="position: absolute; font-family: Arial, Helvetica, sans-serif; font-size: 1.5em; top: 8%; left: 50%; transform:translate(-50%,-50%)">Select <i>FastQ</i> file to upload:</div> 
    <input type="file" name="file" id="file" class="inputUploadfile" accept=".txt, .fq, .FastQ, .fastq" @change="changeFileInfo"/>
    <label for="file"><div style="position: absolute; top: 50%; left: 50%; transform: translate(-50%,-50%);">CLICK TO UPLOAD</div></label>
    <div class="uploadFilenameDisplay">File selected : {{fileName || "None"}} Size : {{fileSize}}</div>

    <input maxlength="100" type="text" id="name" class="box" v-model="labelValue" style="position: absolute; top: 75%; left: 50%; transform: translate(-50%,-50%);" v-if="validLabel" placeholder="Input sequece label / identifier (please don't sql inject)" @change="checkValidLabel">        
    <div v-else style="color: red; font-size: 0.9em; position: absolute; top: 75%; left: 50%; transform: translate(-50%,-50%);" >
        <input maxlength="100" type="password" id="password" class="youDidSomethingBadBox" v-model="formData.password" required placeholder="Password" @change="checkValidLabel"> 
        Password should be atleast 6 characters
    </div>
    <div><button class="uploadSubmitButton" @click="uploadFile">Submit</button></div>
  </div>
  
</template>

<style>


a:link {
  color: white;
  text-decoration: none;
}

a:visited {
  color: white;
  text-decoration: none;
}

a:hover {
  color: white;
}

a:active {
  color: white;
  text-decoration: none;
}

.background{
  position: fixed;
  top: 0px;
  left: 0px;
  padding: 100%;
  background-color: #EBEBEB;
}

.titleWrapper{
    position: fixed;
    top: 0px;
    left: 0px;
    background: #FF6700 ;
    width: 100%;
    padding-top: 1%;
    padding-bottom: 1%;
}

.titleText{
    display: inline;
    padding-left: 1%;
    color: white;
    font-weight: 800;
}

.accountName {
  position: inherit;
  font-family: Arial, Helvetica, sans-serif;
  color: white;
  font-size: 1.5em;
  top: 2%;
  left: 80%;
}

.headingWrapper{
    position: absolute;
    top: 10%;
    left: 50%;
}

.heading{
    transform: translate(-50%);
    font-size: xx-large;
}

.uploadButtonWrapper{
    position: absolute;
    top: 16%;
    left: 50%;
    height: 40%;
    width: 50%;
    transform: translate(-50%);
}

.inputUploadfile {
	width: 0.1px;
	height: 0.1px;
	opacity: 0;
	overflow: hidden;
	position: absolute;
	z-index: -1;
}

.inputUploadfile + label {
    position: absolute;
    left: 50%;
    top: 40%;
    width: 80%;
    height: 50%;
    transform: translate(-50%,-50%);
    text-align: center;
    font-size: 1.7em;
    font-weight: 700;
    color: #EBEBEB;
    background-color: #6a66a3;
    display: inline-block;
    border-radius: 4px;
    border-style: none;
    
}

.inputUploadfile + label:hover{
  cursor: pointer;
}

.uploadSubmitButton {
  position: absolute;
  top: 90%;
  left: 50%;
  transform: translate(-50%);
  outline: none;
  cursor: pointer;
  font-weight: 600;
  border-radius: 3px;
  padding: 12px 24px;
  border: 0;
  color: #3a4149;
  background: #FF6700;
  line-height: 1.15;
  font-size: 1em;

}

.uploadSubmitButton:hover {
  transition: all .1s ease;
  box-shadow: 0 0 0 0 #fff, 0 0 0 3px #1de9b6;
}
.uploadFilenameDisplay{
    position: absolute;
    left: 50%;
    top: 82%;
    transform: translate(-50%);
    white-space: nowrap;
}

.uploadedFileViewer{
    background-color: rgb(161, 61, 99);
    color: #fff;
    position: absolute;
    padding: 2%;
    border-radius: 4px;
    border-style: none;
    white-space: nowrap;
    top: 15%;
    left: 50%;
    transform: translate(-50%);
    font-size: 2em;
}
</style>
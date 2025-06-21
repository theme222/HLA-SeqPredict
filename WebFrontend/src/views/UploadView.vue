<script setup>
/* eslint-disable */
import { ref } from 'vue'
import { reactive } from 'vue'
import axios from 'axios'
import Cookies from 'js-cookie'
import JSZip from "jszip";
import { useRouter, useRoute } from 'vue-router'
import { backendLink } from '@/scripts/AccountFunc'
import { ValidLabel } from '@/scripts/OtherFunc'

const UPLOAD_CHUNK_SIZE = 10 * 1024 * 1024; // 10 MB

const router = useRouter()
const route = useRoute()

const fileName1 = ref(null)
const fileSize1 = ref("0 Bytes")
const fileName2 = ref(null)
const fileSize2 = ref("0 Bytes")
const validLabel = ref(true)
const labelValue = ref('')
const progressUploadedChunks = ref(-1); // This is the value that the progress will look which our interval function will modify based on uploadedChunks
const uploadedChunks = ref(-1); // -1 is our null value
const totalChunks = ref(-1);
let isUploadingFile = false;

const singleEndFile = ref(`${backendLink}/reference/example_files/single_end_example.zip`);
const pairedEndFile = ref(`${backendLink}/reference/example_files/paired_end_example.zip`);

function changeFileInfo(event){
  if (event.target.id == "file1")
  {
    fileName1.value = event.target.files[0].name
    fileSize1.value = formatFileSize(event.target.files[0].size)
  }
  else
  {
    fileName2.value = event.target.files[0].name
    fileSize2.value = formatFileSize(event.target.files[0].size)
  }
}

function formatFileSize(bytes)
{
  const units = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
  let index = 0;
  while (bytes >= 1024 && index < 4) {
    bytes /= 1024;
    index++;
  }
  return `${bytes.toFixed(0)} ${units[index]}`
};

async function uploadFile()
{
  let fileToUpload1 = document.getElementById('file1').files[0]
  let fileToUpload2 = document.getElementById("file2").files[0]
  let fileToUpload;
  let isPaired = false;

  if (!validLabel.value) {return 0}
  if (isUploadingFile) return; 

  if (fileToUpload1 === undefined && fileToUpload2 === undefined){
    alert("No file specified") 
    return 0
  }
  else if (fileToUpload1 && fileToUpload2)
  {
    fileToUpload = new JSZip();
    // add file to zip thingy idk man
    fileToUpload.file("file1.fq", fileToUpload1);
    fileToUpload.file("file2.fq", fileToUpload2);
    fileToUpload = await fileToUpload.generateAsync({type: "blob"});  // Creates zip file
    isPaired = true;
  }
  else if (fileToUpload1) fileToUpload = fileToUpload1;
  else if (fileToUpload2) fileToUpload = fileToUpload2;
  
  isUploadingFile = true;
  uploadedChunks.value = 0;
  totalChunks.value = Math.ceil(fileToUpload.size / UPLOAD_CHUNK_SIZE);

  console.log(isPaired);
  const chunkArray = [];
  for (let i = 0; i < totalChunks.value; i++)
  {
    const formData = new FormData()
    let filenameToSet = labelValue.value || String(Math.floor(Date.now()/1000))
    const chunk = fileToUpload.slice(i * UPLOAD_CHUNK_SIZE, (i+1) * UPLOAD_CHUNK_SIZE);
    formData.append("chunk", chunk, filenameToSet);
    formData.append("chunk_index", i);
    formData.append("total_chunk", totalChunks.value);
    formData.append("is_paired", isPaired);
    chunkArray.push(formData);
  }
  const uploadNextChunk = (index) => {
    if (index >= chunkArray.length) {
      uploadedChunks.value += 1;
      console.log("File upload success");
      setTimeout(()=> router.push("dashboard"), 100);
      return;
    }
    const formData = chunkArray[index];
    axios.post(backendLink+"/api/uploadFile", formData, 
      {headers: {"Content-Type" : 'multipart/form-data', 
                'Authorization': `Bearer ${Cookies.get('session_cookie')}`},
        maxBodyLength: Infinity, // Allow unlimited body size
        maxContentLength: Infinity, // Allow unlimited content length) 
              }
    ).then((response) => {
      console.log(response.data)
      uploadedChunks.value += 1;
      uploadNextChunk(index+1);
    })
    .catch((err) => console.error(err));
  }
  uploadNextChunk(0);

}

function DoProgressBar()
{
  if (uploadedChunks.value == -1) return;
  let alpha = 0.3;
  progressUploadedChunks.value += (uploadedChunks.value - progressUploadedChunks.value) * alpha; 
}

setInterval(DoProgressBar, 50)
</script>

<template>
  <dialog id="infoModal" class="modal absolute">
    <div class="modal-box">
      <h3 class="font-bold text-lg">Hello!</h3>
      <p class="py-4">Your file will be uploaded to our backend servers. Please make sure to read the <RouterLink to="pdpa" class="link">PDPA</RouterLink> for information on how we will use your data. You may remove any data from our servers in the <RouterLink to="profile" class="link"> profile</RouterLink> tab.</p>
      <div class="modal-action">
        <form method="dialog">
          <button class="btn">Close</button>
        </form>
      </div>
    </div>
  </dialog>

  <div class="absolute top-20 grid grid-cols-1 w-full py-10">
    <div class="flex justify-center items-center w-full h-60 md:h-75">
      <div class="w-11/12 max-w-200 h-full shadow-md bg-white rounded-md">

        <div class="flex justify-center items-center w-full h-16 py-12 md:py-18">
          <h1 class="font-bold text-3xl sm:text-4xl">Select file for upload</h1>
        </div>

        <div class="w-full flex justify-center items-center h-20">
          <div class="w-47/50 flex justify-between h-full">
            <div class="fieldset w-24/50">
                <input id="file1" type="file" class="file-input file-input-sm md:file-input-md file-input-bordered file-input-primary w-full" accept=".fq, .FastQ, .fastq .FQ" @change="changeFileInfo"/>
                <span class="fieldset-label">Size : {{fileSize1}}</span>
              </div>
              
              <div class="fieldset w-24/50">
                <input id="file2" type="file" class="file-input file-input-sm md:file-input-md file-input-bordered file-input-secondary w-full" accept=".fq, .FastQ, .fastq .FQ" @change="changeFileInfo"/>
                <span class="fieldset-label">Size : {{fileSize2}}</span>
              </div>
        </div>
      </div>
      


      <div class="w-full flex justify-center">
        <div class="w-47/50 h-10 md:h-12 top-20 flex gap-6 justify-evenly items-center">

          <input v-if="validLabel" v-model="labelValue" @change="validLabel = ValidLabel(labelValue)" type="text" placeholder="Sequence Label Here" class="input input-bordered input-secondary w-9/12 h-full " />
          <input v-else v-model="labelValue" @change="validLabel = ValidLabel(labelValue)" type="text"  placeholder="Sequence Label Here" class="input input-bordered input-error w-9/12 h-full " />

          <div class="h-full flex justify-end items-end"> 
            <button class="btn btn-info btn-outline h-full rounded-md text-sm sm:text-md md:text-xl" onclick="infoModal.showModal()">Info</button>
          </div>

          <button class="btn btn-accent h-full w-32" @click="uploadFile">
            Upload
            <svg class="h-6 w-6" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M17 17H17.01M15.6 14H18C18.9319 14 19.3978 14 19.7654 14.1522C20.2554 14.3552 20.6448 14.7446 20.8478 15.2346C21 15.6022 21 16.0681 21 17C21 17.9319 21 18.3978 20.8478 18.7654C20.6448 19.2554 20.2554 19.6448 19.7654 19.8478C19.3978 20 18.9319 20 18 20H6C5.06812 20 4.60218 20 4.23463 19.8478C3.74458 19.6448 3.35523 19.2554 3.15224 18.7654C3 18.3978 3 17.9319 3 17C3 16.0681 3 15.6022 3.15224 15.2346C3.35523 14.7446 3.74458 14.3552 4.23463 14.1522C4.60218 14 5.06812 14 6 14H8.4M12 15V4M12 4L15 7M12 4L9 7" stroke="#000000" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></svg>
          </button>
        </div>
      </div>
    </div>
  </div>

    <div class="w-full h-30 flex justify-center items-center">
      <div class="w-11/12 max-w-200 bg-white h-16 sm:h-20 rounded-md shadow-md flex justify-center">
        <div class="w-46/50 flex items-center h-full gap-6">
          <h1 class="text-center sm:text-2xl font-semibold">Example Files : </h1>
          <a class="btn btn-soft btn-primary btn-sm md:btn-md" :href="pairedEndFile">Paired End</a>
          <a class="btn btn-soft btn-primary btn-sm md:btn-md" :href="singleEndFile">Single End</a>
        </div> 
      </div>
    </div>
  
    <div v-if="progressUploadedChunks!=-1" class="w-full h-8 flex justify-center items-center">
      <progress :value="progressUploadedChunks" :max="totalChunks" class="progress w-11/12 max-w-200 h-8 transition-all ease-in duration-100 progress-primary"></progress>
    </div>

    </div>


  
</template>
   

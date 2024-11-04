<script setup>
/* eslint-disable */
import { ref } from 'vue'
import { reactive } from 'vue'
import axios from 'axios'
import Cookies from 'js-cookie'
import { useRouter, useRoute } from 'vue-router'
import { backendLink } from '@/scripts/AccountFunc'

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
    const regex = /^[a-zA-Z0-9_]+$/;
    validLabel.value = regex.test(labelValue.value)
}

function uploadFile(){
  console.log(document.getElementById('file').files[0])
  if (!validLabel.value) {return 0}
    if (document.getElementById('file').files[0] === undefined){
      alert("No file specified") 
      return 0
    }
    const formData = new FormData()
    let filenameToSet  = String(Math.floor(Date.now()/1000))
    if (labelValue.value !== "") {filenameToSet = labelValue.value}
    formData.append('file', document.getElementById('file').files[0], filenameToSet)
    axios.post(backendLink+"/api/uploadFile", 
              formData, 
              {headers: {"Content-Type" : 'multipart/form-data', 
                        'Authorization': `Bearer ${Cookies.get('session_cookie')}`} })
    .then(response => {
      alert(response.data)
      router.push({name:'dashboard'})
    })
    .catch(error => {
      alert("Please use alpha numeric characters and _ only")
      console.error('Error uploading file  :', error);
    });
  }

</script>

<template>

  <div class="absolute w-full h-80 top-32 flex justify-center items-center ">
    <div class="w-200 h-full shadow-md bg-white rounded-md">

      <div class="flex justify-center items-center w-full h-16 relative top-8 ">
        <h1 class="font-bold text-4xl">Select file for upload</h1>
      </div>

      <div class="w-full flex justify-evenly relative top-20 h-20">

          <div class="w-8/12">
            <input id="file" type="file" class="file-input file-input-bordered file-input-primary w-full min-w-xs " accept=".txt, .fq, .FastQ, .fastq" @change="changeFileInfo"/>
            <div class="label ">
              <span class="label-text-alt">Size : {{fileSize}}</span>
            </div>
          </div>

          <div class="w-1/4 h-3/5 flex justify-center items-center"> 
            <button class="border border-secondary bg-base-100 font-bold w-full h-full" onclick="infoModal.showModal()">Click me for info</button>
            <dialog id="infoModal" class="modal">
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
          </div>

      </div>
      


      <div class="w-full h-12 relative top-20 flex justify-evenly items-center">

        <input v-if="validLabel" v-model="labelValue" @change="checkValidLabel" type="text" placeholder="File label here (Use alphanumeric and _ only)" class="input input-bordered input-secondary w-9/12 h-full " />
        <input v-else v-model="labelValue" @change="checkValidLabel" type="text"  placeholder="File label here (Use alphanumeric and _ only)" class="input input-bordered input-error w-9/12 h-full " />

        <button class="btn btn-accent h-full w-32" @click="uploadFile">
          Upload
          <svg class="h-6 w-6" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M17 17H17.01M15.6 14H18C18.9319 14 19.3978 14 19.7654 14.1522C20.2554 14.3552 20.6448 14.7446 20.8478 15.2346C21 15.6022 21 16.0681 21 17C21 17.9319 21 18.3978 20.8478 18.7654C20.6448 19.2554 20.2554 19.6448 19.7654 19.8478C19.3978 20 18.9319 20 18 20H6C5.06812 20 4.60218 20 4.23463 19.8478C3.74458 19.6448 3.35523 19.2554 3.15224 18.7654C3 18.3978 3 17.9319 3 17C3 16.0681 3 15.6022 3.15224 15.2346C3.35523 14.7446 3.74458 14.3552 4.23463 14.1522C4.60218 14 5.06812 14 6 14H8.4M12 15V4M12 4L15 7M12 4L9 7" stroke="#000000" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></svg>
        </button>
      </div>

    </div>
  </div>
  
  
</template>


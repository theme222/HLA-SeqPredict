<script setup>
/* eslint-disable */
let func = "This is so that autocomplete doesn't do some bullshit";
import { ref } from "vue";
import { reactive } from "vue";
import axios from "axios";
import Cookies from "js-cookie";
import { useRouter, useRoute } from "vue-router";
import { getAccountInfo, getSequences, backendLink } from "@/scripts/AccountFunc";

const router = useRouter();
const route = useRoute();

const sequenceDropdown = ref([]);
const currentSequence = ref("");
const accountName = ref("");
const accountEmail = ref("");
const deleteButtonStatus = ref("Delete Sequence");

getAccountInfo().then((data) => {
  if (data)
  {
    accountName.value = data["name"];
    accountEmail.value = data["email"];
  }
  else
  {
    router.push({name:"login"})
    alert("Please login first")
  }
});

async function deleteSequence(){
  if (currentSequence.value == "") {return}
  if (deleteButtonStatus.value == "Delete Sequence"){
    deleteButtonStatus.value = "Confirm delete?"
    setTimeout(() => {deleteButtonStatus.value = "Delete Sequence"},3000)
  } else {
    axios.post(backendLink+"/api/deleteSequence", {session_cookie: Cookies.get('session_cookie'), label:currentSequence.value})
    .then(response =>{
      alert("All releted information deleted")
      console.log(response.data)
      sequenceDropdown.value = []
      addToSequenceList()
    })
    .catch(error => {
        console.error('Error getting info :', error);
        alert("Backend server unavailable")

    });
  }
}

async function addToSequenceList() {
  let serverData = await getSequences();
  for (let sequence of serverData) {
    sequenceDropdown.value.push({ value: sequence[1] });
  }
  console.log(sequenceDropdown.value);
}

console.log(accountEmail);
addToSequenceList();
</script>

<template>


  <div class="absolute w-full h-96 top-32 flex justify-center items-center ">
    <div class="w-256 h-full shadow-md bg-white rounded-md">
      
      <div class="relative left-8 w-80 h-24 flex items-center">
        <p class="text-5xl font-bold"> Account info </p>
      </div>

      <div class="w-128 h-64 grid grid-cols-1 grid-rows-3 relative left-8 top-0">

        <label class="form-control w-full">
          <div class="label">
            <span class="label-text">Name</span>
          </div>
          <input  v-model="accountName" @change="checkValidLabel" type="text"  class="input input-bordered input-primary w-full h-16" />
        </label>

        <label class="form-control w-full">
          <div class="label">
            <span class="label-text">Email</span>
          </div>
          <input  v-model="accountEmail" @change="checkValidLabel" type="text" class="input input-bordered input-primary w-full h-16" />
        </label>
        
        <label class="form-control w-full">
          <div class="label">
            <span class="label-text">Password</span>
          </div>
          <input @change="checkValidLabel" type="password" value="idk your password bro it's hashed" class="input input-bordered input-primary w-full h-16" />
        </label>

      </div>

      <div class="w-full h-12 relative -top-12 flex justify-end items-center pointer-events-none">
          <button disabled class="btn btn-accent h-full w-64 relative right-20 pointer-events-auto" @click="uploadFile">
            <p class="px-2">Save changes</p>
            <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-pencil-fill" viewBox="0 0 16 16">
              <path d="M12.854.146a.5.5 0 0 0-.707 0L10.5 1.793 14.207 5.5l1.647-1.646a.5.5 0 0 0 0-.708zm.646 6.061L9.793 2.5 3.293 9H3.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.207zm-7.468 7.468A.5.5 0 0 1 6 13.5V13h-.5a.5.5 0 0 1-.5-.5V12h-.5a.5.5 0 0 1-.5-.5V11h-.5a.5.5 0 0 1-.5-.5V10h-.5a.5.5 0 0 1-.175-.032l-.179.178a.5.5 0 0 0-.11.168l-2 5a.5.5 0 0 0 .65.65l5-2a.5.5 0 0 0 .168-.11z"/>
            </svg> 
          </button>
      </div>

      <div class="w-full relative -top-96 flex justify-end items-center pointer-events-none">
        <div class="avatar placeholder right-20">
          <div class="bg-neutral text-neutral-content rounded-full w-64">
            <span class="text-7xl">{{accountName.toUpperCase().split(' ').map(word => word[0]).join('')}}</span>
          </div>
        </div> 
      </div>

    </div>
  </div>
  
  <div class="absolute w-full h-96 top-144 flex justify-center items-center ">
    <div class="w-256 h-full shadow-md bg-white rounded-md">
      
      <div class="relative left-8 w-80 h-24 flex items-center">
        <p class="text-5xl font-bold"> Sequences </p>
      </div>

      <div class="w-64 h-14 flex justify-start relative left-8 top-0">
        <select v-model="currentSequence" class="select select-primary w-full max-w-xs">
          <option value="" disabled selected>Select a sequence</option>
          <option v-for="i in sequenceDropdown" :key="i" :value="i.value">
            {{ i.value }}
          </option>
        </select>


      </div>

      <div class="w-full h-12 relative left-8 top-40 flex items-center pointer-events-none ">
          <button class="btn btn-error h-full w-64 relative  pointer-events-auto" @click="deleteSequence">
            <p class="px-2">{{deleteButtonStatus}}</p>
            <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-pencil-fill" viewBox="0 0 16 16">
              <path d="M12.854.146a.5.5 0 0 0-.707 0L10.5 1.793 14.207 5.5l1.647-1.646a.5.5 0 0 0 0-.708zm.646 6.061L9.793 2.5 3.293 9H3.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.207zm-7.468 7.468A.5.5 0 0 1 6 13.5V13h-.5a.5.5 0 0 1-.5-.5V12h-.5a.5.5 0 0 1-.5-.5V11h-.5a.5.5 0 0 1-.5-.5V10h-.5a.5.5 0 0 1-.175-.032l-.179.178a.5.5 0 0 0-.11.168l-2 5a.5.5 0 0 0 .65.65l5-2a.5.5 0 0 0 .168-.11z"/>
            </svg> 
          </button>
      </div>

    </div>
  </div>
  
</template>

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
const accountName = ref();
const accountEmail = ref();
const deleteButtonStatus = ref("Delete");

getAccountInfo().then((data) => {
  accountName.value = data["name"];
  accountEmail.value = data["email"];
});

async function deleteSequence(){
  if (deleteButtonStatus.value == "Delete"){
    deleteButtonStatus.value = "Confirm?"
    setTimeout(() => {deleteButtonStatus.value = "Delete"},3000)
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
  <h1 class="heading">Welcome, {{ accountName }} ({{ accountEmail }})</h1>

  <div class="profileSelectDiv">
    <h2 class="heading2">Current sequence :</h2>

    <div
      class="dropdownDiv"
      style="top: 50%; left: 390px; transform: translate(-50%, -50%)"
    >
      <select
        name="sequences"
        id="sequences"
        class="dropdown"
        v-model="currentSequence"
      >
        <option value="" disabled selected>Select a sequence</option>
        <option v-for="i in sequenceDropdown" :key="i" :value="i.value">
          {{ i.value }}
        </option>
      </select>
      <img
        src="@/assets/down-arrow.svg"
        alt="down-arrow"
        style="
          position: absolute;
          filter: invert(99%) sepia(1%) saturate(594%) hue-rotate(271deg)
            brightness(114%) contrast(100%);
          left: 86%;
          top: 50%;
          width: 25%;
          height: 100%;
          pointer-events: none;
          transform: translate(-50%, -50%);
        "
      />
    </div>
    <button class="deleteButton" @click="deleteSequence">
      <img
        src="@/assets/trash-can.svg"
        alt="trash-can"
        style="
          position: absolute;
          height: 45%;
          width: 34%;
          left: 19%;
          top: 50%;
          transform: translate(-50%, -50%);
        "
      />

      <span
        style="
          left: 60%;
          top: 50%;
          transform: translate(-50%, -50%);
          position: absolute;
        "
        >{{deleteButtonStatus}}</span
      >
    </button>
  </div>
</template>

<style>
.profileSelectDiv {
  position: absolute;
  left: 50%;
  top: 20%;
  transform: translate(-50%, -50%);
  width: 100%;
  height: 9%;
}

.heading2 {
  position: absolute;
  top: 35%;
  left: 130px;
  transform: translate(-50%, -50%);
  font-family: Arial, Helvetica, sans-serif;
}

.deleteButton {
  position: absolute;
  top: 50%;
  left: 600px;
  height: 51px;
  width: 140px;
  display: inline-block;
  outline: none;
  cursor: pointer;
  font-weight: 600;
  border-radius: 3px;
  padding: 12px 24px;
  border: 0;
  color: #fff;
  background: #36382e;
  line-height: 1.15;
  font-size: 1em;
  transform: translate(-50%, -50%);
}
.deleteButton:hover {
  background: #B20D30;
  transition: all 0.5s ease;
}
</style>
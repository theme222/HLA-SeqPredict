<script setup>
/* eslint-disable */
let func = "This is so that autocomplete doesn't do some bullshit";
import { ref } from "vue";
import { reactive } from "vue";
import axios from "axios";
import Cookies from "js-cookie";
import { useRouter, useRoute } from "vue-router";
import { getAccountInfo, getSequences, backendLink,Logout } from "@/scripts/AccountFunc";
import { ValidName, ValidEmail, ValidPassword } from "@/scripts/OtherFunc";
import { element } from "plotly.js-dist";
import { sha256 } from "js-sha256";

const router = useRouter();
const route = useRoute();

const sequenceDropdown = ref([]);
const currentSequence = ref("");
const deleteButtonStatus = ref("Delete Sequence");
const changeButtonStatus = ref("Save Changes");

const accountName = ref("");
const accountEmail = ref("");
const accountPassword = ref("");

// default (before changing) value
let _accountName;
let _accountEmail;
let _accountPassword = "";

const nameInput = ref();
const emailInput = ref();
const passwordInput = ref();

// value of element.disabled
const _nameInput = ref(true);
const _emailInput = ref(true);
const _passwordInput = ref(true);


function CheckAll()
{
  let a = 0
  const InvertInput = (element, boolValue) =>
  {
    if (!boolValue)
    {
      a += 1;
      element.classList.remove('input-primary');
      element.classList.add('input-error');
    }
    else
    {
      element.classList.add('input-primary');
      element.classList.remove('input-error');
    }
  }

  if (!_nameInput.value) InvertInput(nameInput.value, ValidName(accountName.value));
  if (!_emailInput.value) InvertInput(emailInput.value, ValidEmail(accountEmail.value));
  if (!_passwordInput.value) InvertInput(passwordInput.value, ValidPassword(accountPassword.value));
  return a == 0;
}

getAccountInfo().then((data) => {
  if (data)
  {
    accountName.value = data["name"];
    accountEmail.value = data["email"];
    _accountName = data["name"];
    _accountEmail = data["email"];
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
    deleteButtonStatus.value = "Confirm Delete?"
    setTimeout(() => {deleteButtonStatus.value = "Delete Sequence"},3000)
  } else {
    axios.post(backendLink+"/api/deleteSequence", {cookie: Cookies.get('session_cookie'), label:currentSequence.value})
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
  if (!serverData) return;
  for (let sequence of serverData) {
    sequenceDropdown.value.push({ value: sequence[1] });
  }
  console.log(sequenceDropdown.value);
}

function ToggleInput(inputName)
{
  switch (inputName){
    case "name":
      _nameInput.value = !nameInput.value.disabled;
      if (!_nameInput.value) accountName.value = _accountName;
      nameInput.value.disabled = !nameInput.value.disabled;
      break;
    case "email":
      _emailInput.value = !emailInput.value.disabled;
      if (!_emailInput.value) accountEmail.value = _accountEmail;
      emailInput.value.disabled = !emailInput.value.disabled;
      break;
    case "password":
      _passwordInput.value = !passwordInput.value.disabled;
      if (!_passwordInput.value) accountPassword.value = _accountPassword;
      passwordInput.value.disabled = !passwordInput.value.disabled;
      break;
  }
}

async function ChangeAccountInfo()
{
  if (_nameInput.value && _emailInput.value && _passwordInput.value) return;
  if (!CheckAll()) return;
  if (changeButtonStatus.value == "Save Changes")
  {
    changeButtonStatus.value = "Are you sure?";
    setTimeout(() => {changeButtonStatus.value = "Save Changes"},3000)
    return;
  }
  
  let dict = {};
  if (!_nameInput.value) dict["name"] = accountName.value;
  if (!_emailInput.value) dict["email"] = accountEmail.value;
  if (!_passwordInput.value) dict["password"] = sha256(accountPassword.value);
  dict["cookie"] = Cookies.get('session_cookie');

  axios.post(backendLink+"/api/changeAccountInfo", dict)
  .then(response =>{
    alert("Saved changes");
    window.location.reload();
  })
  .catch(error => {
    if (error.response['data'] == "Duplicate email") alert("Duplicate email");
    else
    {
      console.error('Error getting info :', error);
      alert("Backend server unavailable")
    }
  });


}

console.log(accountEmail);
addToSequenceList();
</script>

<template>


  <div class="absolute w-full h-96 top-32 flex justify-center items-center ">
    <div class="w-256 h-full shadow-md bg-white rounded-md">
      
      <div class="relative left-8 w-128 h-24 flex items-center justify-between ">
        <p class="text-5xl font-bold"> Account info </p>
        <div class="btn w-40 h-14 text-xl" @click="Logout();router.push({name: 'login'})">Log Out 
          <svg class="w-5 h-6" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" width="256" height="256" viewBox="0 0 256 256" xml:space="preserve">
          <defs>
          </defs>
          <g style="stroke: none; stroke-width: 0; stroke-dasharray: none; stroke-linecap: butt; stroke-linejoin: miter; stroke-miterlimit: 10; fill: none; fill-rule: nonzero; opacity: 1;" transform="translate(1.4065934065934016 1.4065934065934016) scale(2.81 2.81)" >
            <path d="M 86.356 46.27 c 0.031 -0.065 0.059 -0.131 0.085 -0.199 c 0.042 -0.11 0.076 -0.222 0.104 -0.336 c 0.016 -0.062 0.034 -0.123 0.046 -0.186 c 0.034 -0.181 0.055 -0.364 0.055 -0.548 l 0 0 c 0 0 0 0 0 0 c 0 -0.184 -0.022 -0.367 -0.055 -0.548 c -0.012 -0.064 -0.03 -0.124 -0.046 -0.186 c -0.029 -0.114 -0.062 -0.226 -0.104 -0.336 c -0.026 -0.068 -0.055 -0.134 -0.086 -0.199 c -0.046 -0.099 -0.099 -0.194 -0.156 -0.288 c -0.039 -0.063 -0.077 -0.126 -0.12 -0.186 c -0.02 -0.027 -0.033 -0.057 -0.054 -0.084 L 74.316 27.93 c -1.009 -1.313 -2.894 -1.561 -4.207 -0.551 c -1.313 1.009 -1.561 2.893 -0.551 4.207 L 77.56 42 H 30.903 c -1.657 0 -3 1.343 -3 3 c 0 1.657 1.343 3 3 3 h 46.656 l -8.001 10.414 c -1.01 1.314 -0.763 3.197 0.551 4.207 c 0.545 0.419 1.188 0.621 1.826 0.621 c 0.9 0 1.79 -0.403 2.381 -1.172 l 11.71 -15.242 c 0.021 -0.027 0.035 -0.057 0.055 -0.085 c 0.043 -0.06 0.08 -0.122 0.119 -0.184 C 86.257 46.464 86.31 46.369 86.356 46.27 z" style="stroke: none; stroke-width: 1; stroke-dasharray: none; stroke-linecap: butt; stroke-linejoin: miter; stroke-miterlimit: 10; fill: rgb(0,0,0); fill-rule: nonzero; opacity: 1;" transform=" matrix(1 0 0 1 0 0) " stroke-linecap="round" />
            <path d="M 60.442 90 H 9.353 c -1.657 0 -3 -1.343 -3 -3 V 3 c 0 -1.657 1.343 -3 3 -3 h 51.089 c 1.657 0 3 1.343 3 3 v 30.054 c 0 1.657 -1.343 3 -3 3 s -3 -1.343 -3 -3 V 6 H 12.353 v 78 h 45.089 V 55.61 c 0 -1.657 1.343 -3 3 -3 s 3 1.343 3 3 V 87 C 63.442 88.657 62.1 90 60.442 90 z" style="stroke: none; stroke-width: 1; stroke-dasharray: none; stroke-linecap: butt; stroke-linejoin: miter; stroke-miterlimit: 10; fill: rgb(0,0,0); fill-rule: nonzero; opacity: 1;" transform=" matrix(1 0 0 1 0 0) " stroke-linecap="round" />
          </g>
          </svg>
        </div>
      </div>

      <div class="w-128 h-64 grid grid-cols-1 grid-rows-3 relative left-8 top-0">

        <!-- Name input -->
        <label class="form-control w-full">
          <div class="label">
            <span class="label-text">Name</span>
          </div>
          <div class="w-full h-16 flex justify-center items-center gap-2">
            <input ref="nameInput" id="nameInput" v-model="accountName" disabled  type="text" @change="CheckAll" class="input input-bordered input-primary w-full h-full" />

            <div v-if="_nameInput" class="tooltip" data-tip="Edit">
              <div @click="ToggleInput('name')" class="btn btn-outline btn-warning w-12">            
                  <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-pencil-fill" viewBox="0 0 16 16">
                  <path d="M12.854.146a.5.5 0 0 0-.707 0L10.5 1.793 14.207 5.5l1.647-1.646a.5.5 0 0 0 0-.708zm.646 6.061L9.793 2.5 3.293 9H3.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.207zm-7.468 7.468A.5.5 0 0 1 6 13.5V13h-.5a.5.5 0 0 1-.5-.5V12h-.5a.5.5 0 0 1-.5-.5V11h-.5a.5.5 0 0 1-.5-.5V10h-.5a.5.5 0 0 1-.175-.032l-.179.178a.5.5 0 0 0-.11.168l-2 5a.5.5 0 0 0 .65.65l5-2a.5.5 0 0 0 .168-.11z"/>
                </svg> 
              </div>
            </div>            
            <div v-else class="tooltip" data-tip="Cancel">
              <div @click="ToggleInput('name')" class="btn btn-outline btn-error w-12">            
                    <svg version="1.0" xmlns="http://www.w3.org/2000/svg" class="w-6 h-6" 
                    width="1280.000000pt" height="1280.000000pt" viewBox="0 0 1280.000000 1280.000000"
                    preserveAspectRatio="xMidYMid meet">
                    <g transform="translate(0.000000,1280.000000) scale(0.100000,-0.100000)"
                    fill="currentColor" stroke="none">
                    <path d="M2717 12790 c-85 -15 -190 -51 -272 -92 -78 -39 -100 -60 -1128
                    -1092 -576 -578 -1068 -1078 -1093 -1110 -66 -86 -92 -134 -124 -229 -97 -286
                    -43 -606 153 -900 53 -80 245 -277 1477 -1516 778 -783 1415 -1427 1415 -1431
                    0 -4 -631 -634 -1402 -1398 -772 -765 -1428 -1419 -1458 -1453 -324 -363 -376
                    -885 -121 -1213 50 -64 2041 -2024 2146 -2113 166 -140 425 -208 661 -174 179
                    26 383 114 525 227 37 30 701 684 1476 1455 l1410 1402 227 -229 c125 -126
                    504 -508 841 -849 1599 -1615 1798 -1812 1905 -1884 112 -75 229 -128 361
                    -162 74 -20 113 -24 234 -23 127 1 156 4 233 27 109 34 185 71 262 132 89 69
                    2110 2120 2157 2189 220 321 186 767 -86 1126 -41 54 -582 607 -1466 1499
                    l-1402 1414 84 84 c45 46 690 686 1432 1421 891 884 1369 1365 1410 1419 106
                    141 182 314 211 482 40 228 -9 453 -139 635 -42 59 -2053 2064 -2152 2146
                    -129 106 -322 170 -517 170 -224 0 -478 -92 -677 -245 -35 -27 -696 -680
                    -1471 -1452 l-1408 -1404 -653 658 c-359 362 -931 939 -1272 1283 -891 900
                    -940 948 -1038 1013 -187 125 -368 185 -573 192 -66 2 -142 0 -168 -5z"/>
                    </g>
                    </svg>
              </div>
            </div>

          </div>
        </label>
        <!-- Name input -->

        <!-- Email input -->
        <label class="form-control w-full">
          <div class="label">
            <span class="label-text">Email</span>
          </div>
          <div class="w-full h-16 flex justify-center items-center gap-2">
            <input ref="emailInput" id="emailInput" v-model="accountEmail" disabled  type="text" @change="CheckAll" class="input input-bordered input-primary w-full h-full" />

            <div v-if="_emailInput" class="tooltip" data-tip="Edit">
              <div @click="ToggleInput('email')" class="btn btn-outline btn-warning w-12">            
                  <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-pencil-fill" viewBox="0 0 16 16">
                  <path d="M12.854.146a.5.5 0 0 0-.707 0L10.5 1.793 14.207 5.5l1.647-1.646a.5.5 0 0 0 0-.708zm.646 6.061L9.793 2.5 3.293 9H3.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.207zm-7.468 7.468A.5.5 0 0 1 6 13.5V13h-.5a.5.5 0 0 1-.5-.5V12h-.5a.5.5 0 0 1-.5-.5V11h-.5a.5.5 0 0 1-.5-.5V10h-.5a.5.5 0 0 1-.175-.032l-.179.178a.5.5 0 0 0-.11.168l-2 5a.5.5 0 0 0 .65.65l5-2a.5.5 0 0 0 .168-.11z"/>
                </svg> 
              </div>
            </div>            
            <div v-else class="tooltip" data-tip="Cancel">
              <div @click="ToggleInput('email')" class="btn btn-outline btn-error w-12">            
                    <svg version="1.0" xmlns="http://www.w3.org/2000/svg" class="w-6 h-6" 
                    width="1280.000000pt" height="1280.000000pt" viewBox="0 0 1280.000000 1280.000000"
                    preserveAspectRatio="xMidYMid meet">
                    <g transform="translate(0.000000,1280.000000) scale(0.100000,-0.100000)"
                    fill="currentColor" stroke="none">
                    <path d="M2717 12790 c-85 -15 -190 -51 -272 -92 -78 -39 -100 -60 -1128
                    -1092 -576 -578 -1068 -1078 -1093 -1110 -66 -86 -92 -134 -124 -229 -97 -286
                    -43 -606 153 -900 53 -80 245 -277 1477 -1516 778 -783 1415 -1427 1415 -1431
                    0 -4 -631 -634 -1402 -1398 -772 -765 -1428 -1419 -1458 -1453 -324 -363 -376
                    -885 -121 -1213 50 -64 2041 -2024 2146 -2113 166 -140 425 -208 661 -174 179
                    26 383 114 525 227 37 30 701 684 1476 1455 l1410 1402 227 -229 c125 -126
                    504 -508 841 -849 1599 -1615 1798 -1812 1905 -1884 112 -75 229 -128 361
                    -162 74 -20 113 -24 234 -23 127 1 156 4 233 27 109 34 185 71 262 132 89 69
                    2110 2120 2157 2189 220 321 186 767 -86 1126 -41 54 -582 607 -1466 1499
                    l-1402 1414 84 84 c45 46 690 686 1432 1421 891 884 1369 1365 1410 1419 106
                    141 182 314 211 482 40 228 -9 453 -139 635 -42 59 -2053 2064 -2152 2146
                    -129 106 -322 170 -517 170 -224 0 -478 -92 -677 -245 -35 -27 -696 -680
                    -1471 -1452 l-1408 -1404 -653 658 c-359 362 -931 939 -1272 1283 -891 900
                    -940 948 -1038 1013 -187 125 -368 185 -573 192 -66 2 -142 0 -168 -5z"/>
                    </g>
                    </svg>
              </div>
            </div>

          </div>
        </label>
        <!-- Email input -->

        <!-- Password input -->
        <label class="form-control w-full">
          <div class="label">
            <span class="label-text">Password</span>
          </div>
          <div class="w-full h-16 flex justify-center items-center gap-2">
            <input ref="passwordInput" id="passwordInput" v-model="accountPassword" disabled type="password" @change="CheckAll" class="input input-bordered input-primary w-full h-full" />

            <div v-if="_passwordInput" class="tooltip" data-tip="Edit">
              <div @click="ToggleInput('password')" class="btn btn-outline btn-warning w-12">            
                  <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-pencil-fill" viewBox="0 0 16 16">
                  <path d="M12.854.146a.5.5 0 0 0-.707 0L10.5 1.793 14.207 5.5l1.647-1.646a.5.5 0 0 0 0-.708zm.646 6.061L9.793 2.5 3.293 9H3.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.5h.5a.5.5 0 0 1 .5.5v.207zm-7.468 7.468A.5.5 0 0 1 6 13.5V13h-.5a.5.5 0 0 1-.5-.5V12h-.5a.5.5 0 0 1-.5-.5V11h-.5a.5.5 0 0 1-.5-.5V10h-.5a.5.5 0 0 1-.175-.032l-.179.178a.5.5 0 0 0-.11.168l-2 5a.5.5 0 0 0 .65.65l5-2a.5.5 0 0 0 .168-.11z"/>
                </svg> 
              </div>
            </div>            
            <div v-else class="tooltip" data-tip="Cancel">
              <div @click="ToggleInput('password')" class="btn btn-outline btn-error w-12">            
                    <svg version="1.0" xmlns="http://www.w3.org/2000/svg" class="w-6 h-6" 
                    width="1280.000000pt" height="1280.000000pt" viewBox="0 0 1280.000000 1280.000000"
                    preserveAspectRatio="xMidYMid meet">
                    <g transform="translate(0.000000,1280.000000) scale(0.100000,-0.100000)"
                    fill="currentColor" stroke="none">
                    <path d="M2717 12790 c-85 -15 -190 -51 -272 -92 -78 -39 -100 -60 -1128
                    -1092 -576 -578 -1068 -1078 -1093 -1110 -66 -86 -92 -134 -124 -229 -97 -286
                    -43 -606 153 -900 53 -80 245 -277 1477 -1516 778 -783 1415 -1427 1415 -1431
                    0 -4 -631 -634 -1402 -1398 -772 -765 -1428 -1419 -1458 -1453 -324 -363 -376
                    -885 -121 -1213 50 -64 2041 -2024 2146 -2113 166 -140 425 -208 661 -174 179
                    26 383 114 525 227 37 30 701 684 1476 1455 l1410 1402 227 -229 c125 -126
                    504 -508 841 -849 1599 -1615 1798 -1812 1905 -1884 112 -75 229 -128 361
                    -162 74 -20 113 -24 234 -23 127 1 156 4 233 27 109 34 185 71 262 132 89 69
                    2110 2120 2157 2189 220 321 186 767 -86 1126 -41 54 -582 607 -1466 1499
                    l-1402 1414 84 84 c45 46 690 686 1432 1421 891 884 1369 1365 1410 1419 106
                    141 182 314 211 482 40 228 -9 453 -139 635 -42 59 -2053 2064 -2152 2146
                    -129 106 -322 170 -517 170 -224 0 -478 -92 -677 -245 -35 -27 -696 -680
                    -1471 -1452 l-1408 -1404 -653 658 c-359 362 -931 939 -1272 1283 -891 900
                    -940 948 -1038 1013 -187 125 -368 185 -573 192 -66 2 -142 0 -168 -5z"/>
                    </g>
                    </svg>
              </div>
            </div>

          </div>
        </label>
        <!-- Password input -->

      </div>

      <div class="w-full h-12 relative -top-12 flex justify-end items-center pointer-events-none">
          <button class="btn btn-accent h-full w-64 relative right-20 pointer-events-auto" @click="ChangeAccountInfo">
            <p class="px-2">{{ changeButtonStatus }}</p>
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

<script setup>
/* eslint-disable */
import { getAccountInfo, getSequences, backendLink} from "@/scripts/AccountFunc";
import { ref } from "vue";
import axios from "axios";
import router from "@/router";
import { ErrorHandler, ScrollTo } from "@/scripts/OtherFunc";
import { sharedData } from "@/scripts/SharedData";


const sequenceDropdown = ref([]);
const selectedSequence = ref("");

const logModal = ref();
const logTerminal = ref();
const logPages = ref();
const currentPageNumber = ref(1);
const currentTotalLogPage = ref(1);
let currentSequenceLog;
let logIsOpen = false;

const currentSequence = ref("");
const deleteButtonStatus = ref("Delete Sequence");

function gotoProfile()
{
  sharedData["patientName"] = "test";
  router.push({name:"profile"})
}

async function GetReport() 
{
  if (!currentSequence.value) return;
  let seqId = currentSequence.value[0];
  let label = currentSequence.value[1];
  let response = await axios
      .post(backendLink+"/api/getResults/hla_la", {
        cookie: Cookies.get("session_cookie"),
        label: label,
        sequence_id: seqId,
      })
      .catch(console.error);
  sharedData.sequenceLabel = label
  sharedData.sequenceId = seqId
  sharedData.drugList = response.data
  const genotypeList = new Set();
  for (let item of response.data)
  {
    genotypeList.add(item.allele);
  }
  sharedData.genotypeList = Array.from(genotypeList);
  router.push("results")
}

async function RunHLA_LA()
{
  if (!currentSequence.value) return;
  let seqId = currentSequence.value[0];
  let label = currentSequence.value[1];
  let response = await axios
      .post(backendLink+"/api/run/hla_la", {
        cookie: Cookies.get("session_cookie"),
        label: label,
        sequence_id: seqId
      })
      .catch(ErrorHandler);
  if (!response) return;
  alert(response.data)
}

function RunIGV()
{
  if (!currentSequence.value) return;
  router.push({name:"igv", query: {sequence_id: currentSequence.value[0]}});
}

async function GetSequenceList() {
  let serverData = await getSequences();
  if (!serverData) return;
  for (let sequence of serverData) {
    sequenceDropdown.value.push(sequence);
  }
  console.log(sequenceDropdown.value);
}

function ChangeIgvTrackOptions()
{
  window.sharedData.igvShowAlignments.boolVal = !window.sharedData.igvShowAlignments.boolVal;
  window.sharedData.igvShowAlignments.completed = false;
}

async function DeleteSequence(){
  if (!currentSequence.value) return;
  if (deleteButtonStatus.value == "Delete Sequence"){
    deleteButtonStatus.value = "Confirm Delete?"
    setTimeout(() => {deleteButtonStatus.value = "Delete Sequence"},3000)
  } else {
    console.log(currentSequence.value);
    
    axios.post(backendLink+"/api/deleteSequence", {cookie: Cookies.get('session_cookie'), sequence_id:currentSequence.value[0]})
    .then(response =>{
      alert("All releted information deleted")
      console.log(response.data)
      sequenceDropdown.value = []
      GetSequenceList()
    })
    .catch(error => {
        console.error('Error getting info :', error);
        alert("Backend server unavailable")

    }); 
  }
}

function OpenLogViewer()
{
  if (!currentSequence.value) return;
  logIsOpen = true;
  logModal.value.showModal();
  logTerminal.value.innerHTML = "";
  currentSequenceLog = '';
  GetLogs();
}

function DisplayLogPage()
{
  logTerminal.value.innerHTML = "";
  let count = (currentPageNumber.value-1) * 25 + 1;

  let logText = currentSequenceLog.split("\n").slice(25*(currentPageNumber.value-1),25*currentPageNumber.value)
  for (let line of logText)
  {
    if (line[0] == ">") logTerminal.value.innerHTML += `<pre data-prefix='>' class='text-accent'><code>${line.slice(2)}</code></pre>`;
    else if (line[0] == "!") logTerminal.value.innerHTML += `<pre data-prefix='!' class='text-error'><code>${line.slice(2)}</code></pre>`;
    else logTerminal.value.innerHTML += `<pre data-prefix="${count}" class='text-base-300'><code>${line}</code></pre>` ;
    count += 1;
  }
}


async function GetLogs() // I'd like to personally say sorry to all the websocket believers.
{
  if (!logIsOpen) return;
  if (!currentSequence.value) return;
  if (!logTerminal.value || !logModal.value) return;
  let response = await axios.post(backendLink+"/api/getSequenceLog", 
    {
      cookie: Cookies.get("session_cookie"),
      sequence_id: currentSequence.value[0] 
    }
  ).catch(ErrorHandler)
  
  currentSequenceLog = response.data.logs;
  currentTotalLogPage.value = Math.floor(currentSequenceLog.split("\n").length / 25);
  DisplayLogPage();
}

function SetPageNumber(value)
{
  currentPageNumber.value = parseInt(value,10);
  const clamp = (num, min, max) => Math.min(Math.max(num, min), max);
  currentPageNumber.value = clamp(currentPageNumber.value, 1, currentTotalLogPage.value); 
  DisplayLogPage();
}

function CloseLogViewer()
{
  logIsOpen = false;
  currentSequenceLog = '';
  currentTotalLogPage.value = 1;
  currentPageNumber.value = 1;
}


GetSequenceList();
setInterval(GetLogs, 1000);
</script>

<template>

  <dialog class="modal absolute" ref="logModal">
    <div class="modal-box max-w-11/12">
      <form method="dialog">
        <button class="btn btn-sm btn-circle btn-ghost absolute right-2 top-2" @click="CloseLogViewer">✕</button>
      </form>
      <div class="w-full text-center text-3xl font-bold pb-5"> Logs </div>
      <div class="mockup-code w-full" ref="logTerminal">
        <pre data-prefix=">" class="text-primary"><code>ls -lh /app</code></pre>
        <pre><code>installing...</code></pre>
        <pre data-prefix="3" class="bg-warning text-warning-content"><code>Error!</code></pre>
      </div>
      <!-- Probably one of the best code I've ever written in my life-->
      <div class="w-full h-16 flex justify-center items-end gap-4" ref="logPages">
        <div class="btn btn-ghost w-10 h-10 text-center text-xl" @click="SetPageNumber(currentPageNumber-1)">🠈</div>
        <div v-if="currentPageNumber != 1" class="btn btn-ghost w-10 h-10 text-center text-xl" @click="SetPageNumber(1)">1</div>
        <input v-if="currentPageNumber - 1 > 1" v-model="pageInputSelector1" type="text" placeholder="..." class="input input-ghost w-14 h-10 text-center text-sm" @keyup.enter="SetPageNumber(pageInputSelector1)"/>
        <div class="btn btn-primary w-10 h-10 text-center text-xl">{{ currentPageNumber }}</div>
        <input v-if="currentTotalLogPage - currentPageNumber > 1" v-model="pageInputSelector2" type="text" placeholder="..." class="input input-ghost w-14 h-10 text-center text-sm" @keyup.enter="SetPageNumber(pageInputSelector2)"/>  
        <div v-if="currentTotalLogPage != currentPageNumber" class="btn btn-ghost w-10 h-10 text-center text-xl" @click="SetPageNumber(currentTotalLogPage)">{{currentTotalLogPage}}</div>
        <div class="btn btn-ghost w-10 h-10 text-center text-xl" @click="SetPageNumber(currentPageNumber+1)">🠊</div>
      </div>

    </div>
  </dialog>

  <div id="selectedSequence" class="hidden">{{selectedSequence}}</div>

  <div class="absolute top-20 grid gap-8 grid-cols-1 w-full pt-10 pb-40">
  <div class="flex justify-center items-center h-150 md:h-200 w-full">
    <div class="w-11/12 max-w-300 min-w-96 h-full shadow-md bg-white rounded-md">
      <div class="flex items-center justify-center w-full h-8/50 inherit">
        <h1 class="text-4xl font-bold ">Dashboard</h1>
      </div>
      <div class="flex justify-center h-8/10 ">
        <div class="overflow-x-auto w-11/12 ">
          <table class="table table-pin-rows table-sm lg:table-md">
            <thead>
              <tr>
                <th>ID</th>
                <th>Label</th>
                <th class="hidden sm:table-cell">Date Uploaded</th>
                <th class="hidden md:table-cell">Status</th>
                <th>Details</th>
              </tr>
            </thead>
            <tbody>
              
              <tr v-for="i in sequenceDropdown" :key="i" :value="i[0]">
                <th>{{i[0]}}</th>  <!-- Id -->
                <th class="text-primary">{{i[1]}}</th>  <!-- Label -->
                <th class="text-secondary hidden sm:table-cell">{{i[2]}}</th>  <!-- Date -->
                <th class="text-lg hidden md:table-cell">{{i[3]}}</th>  <!-- Status -->
                <th>
                  <button class="btn btn-secondary btn-md rounded-4xl sm:rounded-md" @click="currentSequence = i; ScrollTo('detailsSection')">
                    <span class="px-2 hidden sm:block">Details</span>   
                    <img src="@/assets/down-arrow-white.svg" alt="" class="w-12 h-10 sm:w-8 sm:h-8">
                  </button>
                </th> 
              </tr>

            </tbody>
          </table>
        </div>  
      </div>
      
    </div>  
  </div>

  <div id="detailsSection" class="w-full h-20 flex justify-center">
      <div class="w-11/12 max-w-300 bg-secondary h-full flex items-center justify-center rounded-md">

        <div class="w-48/50 flex items-center justify-between">
          <div class="w-1/4 h-14 flex justify-start">
            <select v-model="currentSequence" class="select select-sm md:select-md select-secondary w-full h-full text-xl">
              <option value="" disabled selected class="">Select a sequence</option>
              <option v-for="i in sequenceDropdown" :key="i[0]" :value="i" class="">
                {{ `${i[0]} : ${i[1]}` }}
              </option>
            </select>
          </div>

          <div class="h-full flex justify-center items-center">
            <div class="indicator">
              <span v-if="currentSequence && currentSequence[3] == 'IDLE'" class="w-4 h-4 animate-ping indicator-item status status-success"></span>
              <span v-if="currentSequence && currentSequence[3] == 'IDLE'" class="w-3 h-3 indicator-item status status-success"></span>
              <span v-if="currentSequence && currentSequence[3] != 'IDLE'" class="w-4 h-4 animate-ping indicator-item status status-warning"></span>
              <span v-if="currentSequence && currentSequence[3] != 'IDLE'" class="w-3 h-3 indicator-item status status-warning"></span>
              <div class="md:text-xl font-semibold flex justify-center items-center bg-base-200 shadow-md h-14 px-10 rounded-sm">
                Status : {{ currentSequence[3] }} 
              </div>
            </div>
          </div>

        <div class="w-1/4 flex justify-end items-center">
          <button class="btn bg-secondary-content w-full h-14" @click="OpenLogViewer">
            <p class="md:px-2 lg:text-lg"><span class="hidden md:inline">View</span> Logs </p>
            <img src="@/assets/terminal.svg" alt="Icon" class="w-7 h-5 md:w-8 md:h-8">
          </button>
        </div>
      </div>


    </div>

  </div>
  <div class="w-full h-100 flex justify-center items-center">
    <div class="w-11/12 max-w-300 h-full rounded-md gap-10 grid grid-cols-1 lg:grid-cols-2">
      <div class="border-primary rounded-sm shadow-lg h-60 border-2 bg-base-200">
        <span class="w-full h-20 flex justify-center items-center font-semibold text-3xl">Info</span>
        <div class="w-full text-neutral *:px-5 text-xl grid grid-cols-1">
          <p>Sequence ID : {{ currentSequence[0] }}</p>
          <p>Sequence Label : {{ currentSequence[1] }}</p>
          <p>Upload Time : {{ currentSequence[2] }}</p>
          <p>Ran for IGV : {{ currentSequence[4] }}</p>
          <p>Ran for HLA*LA : {{ currentSequence[5] }} </p>
        </div>
      </div>

      <div class="border-primary rounded-sm shadow-lg h-60 border-2 bg-base-200">
        <span class="w-full h-20 flex justify-center items-center font-semibold text-3xl">Actions</span>
        <div class="w-full flex justify-center items-center">
          <div class="w-11/12 text-neutral grid grid-cols-2 lg:grid-cols-2 gap-6">
            <button :class="{'btn-disabled': !currentSequence}" class="btn bg-cyan-200" @click="RunIGV">
              <p class="px-2"> Visualize IGV </p>
              <img src="@/assets/eye-icon.svg" alt="Icon" class="w-5 h-5">
            </button>
            <button :class="{'btn-disabled': !currentSequence}" class="btn bg-amber-200" @click="RunHLA_LA">
              <p class="px-2"> Run HLA*LA </p>
              <img src="@/assets/dna.svg" alt="Icon" class="w-7 h-7">
            </button>
            <button :class="{'btn-disabled': !currentSequence}" class="btn bg-purple-300" @click="GetReport">
              <p class="px-2"> Get Report </p>
              <img src="@/assets/report.svg" alt="Icon" class="w-6 h-6">
            </button>
            <button :class="{'btn-disabled': !currentSequence}" class="btn bg-rose-400" @click="DeleteSequence">
              <p class="px-2"> {{ deleteButtonStatus }} </p>
              <img src="@/assets/trash-can.svg" alt="Icon" class="w-5 h-6">
            </button>
          </div>
        </div>
     </div>
 
    </div>
  </div>
</div>
</template>

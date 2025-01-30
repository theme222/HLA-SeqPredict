<script setup>
/* eslint-disable */
import { getAccountInfo, getSequences, backendLink} from "@/scripts/AccountFunc";
import { ref } from "vue";
import axios from "axios";
import router from "@/router";
import { ScrollTo } from "@/scripts/OtherFunc";
import { sharedData } from "@/scripts/SharedData";


const sequenceDropdown = ref([]);
const selectedSequence = ref("");
const typingResults = ref({
  hla_la: ["", "", "", ""],
  optitype: ["", "", "", ""],
  hisat_genotype: ["", "", "", ""],
  snp_bridge: ["", "", "", ""],
});

function gotoProfile()
{
  sharedData["patientName"] = "test";
  router.push({name:"profile"})
}

async function getReport(seq_id,label) 
{
  let response = await axios
      .post(backendLink+"/api/getResults/hla_la", {
        cookie: Cookies.get("session_cookie"),
        label: label,
      })
      .catch(console.error);
  sharedData.sequenceLabel = label
  sharedData.sequenceId = seq_id
  sharedData.drugList = response.data
  const genotypeList = new Set();
  for (let item of response.data)
  {
    genotypeList.add(item.allele);
  }
  sharedData.genotypeList = Array.from(genotypeList);
  router.push("results")
}

async function runTool(sequence_label)
{
  let response = await axios
      .post(backendLink+"/api/run/hla_la", {
        cookie: Cookies.get("session_cookie"),
        label: sequence_label,
      })
      .catch(console.error);
  alert(response.data)
}

async function addToSequenceList() {
  let serverData = await getSequences();
  let labelIndex = 1;
  if (!serverData) return;
  for (let sequence of serverData) {
    sequenceDropdown.value.push(sequence);
  }
  console.log(sequenceDropdown.value);
}




addToSequenceList();
</script>

<template>

  <div id="selectedSequence" class="hidden">{{selectedSequence}}</div>

  <div class="absolute top-20 grid gap-16 grid-cols-1 w-full py-10">
  <div class="flex justify-center items-center h-200 w-full">
    <div class="w-11/12 min-w-96 h-full shadow-md bg-white rounded-md">
      <div class="flex items-center justify-center w-full h-32 inherit">
        <h1 class="text-4xl font-bold ">Dashboard</h1>
      </div>
      <div class="flex justify-center h-4/5 ">
        <div class="overflow-x-auto w-11/12 ">
          <table class="table table-pin-rows  ">
            <thead>
              <tr>
                <th>ID</th>
                <th>Label</th>
                <th>Date Uploaded</th>
                <th>HLA*LA</th>
                <th>View</th>
                <th>Report</th>
                <th>Delete</th>
              </tr>
            </thead>
            <tbody>
              
              <tr v-for="i in sequenceDropdown" :key="i" :value="i[1]">
                <th>{{i[0]}}</th>
                <th>{{i[1]}}</th>
                <th>{{i[2]}}</th>
                
                <th>
                  <button class="btn bg-blue-100 border-blue-100 " @click="runTool(i[1])"><p class="px-2"> Run </p>
                  </button>
                </th>

                <th>
                  <button class="btn bg-black-100 border-black-100" @click="ScrollTo('igvSection'); selectedSequence = i[1]"><p class="px-2"> View </p>
                    <svg fill="#000" class="h-5 w-5" version="1.1" id="Capa_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" 
                      viewBox="0 0 488.85 488.85" xml:space="preserve">
                    <g>
                      <path d="M244.425,98.725c-93.4,0-178.1,51.1-240.6,134.1c-5.1,6.8-5.1,16.3,0,23.1c62.5,83.1,147.2,134.2,240.6,134.2
                        s178.1-51.1,240.6-134.1c5.1-6.8,5.1-16.3,0-23.1C422.525,149.825,337.825,98.725,244.425,98.725z M251.125,347.025
                        c-62,3.9-113.2-47.2-109.3-109.3c3.2-51.2,44.7-92.7,95.9-95.9c62-3.9,113.2,47.2,109.3,109.3
                        C343.725,302.225,302.225,343.725,251.125,347.025z M248.025,299.625c-33.4,2.1-61-25.4-58.8-58.8c1.7-27.6,24.1-49.9,51.7-51.7
                        c33.4-2.1,61,25.4,58.8,58.8C297.925,275.625,275.525,297.925,248.025,299.625z"/>
                    </g>
                    </svg>
                </button>
              </th>

                <th>
                  <button class="btn bg-teal-200 border-teal-200" @click="getReport(i[0],i[1])"><p class="px-2"> Report </p>
                    <svg fill="#000000" class="h-6 w-6" version="1.1" id="Capa_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" 
                      viewBox="0 0 612 612" xml:space="preserve">
                    <g>
                      <g>
                        <g>
                          <g>
                            <path d="M577.661,612H34.339c-10.358,0-18.751-8.396-18.751-18.751V18.751C15.588,8.396,23.98,0,34.339,0h543.322
                              c10.355,0,18.751,8.396,18.751,18.751v574.497C596.412,603.604,588.016,612,577.661,612z M53.09,574.497h505.82V37.502H53.09
                              V574.497z"/>
                          </g>
                          <g>
                            <path d="M476.951,157.596H135.047c-10.355,0-18.751-8.393-18.751-18.751c0-10.355,8.396-18.751,18.751-18.751h341.905
                              c10.355,0,18.751,8.396,18.751,18.751C495.702,149.204,487.307,157.596,476.951,157.596z"/>
                          </g>
                          <g>
                            <path d="M476.951,269.033H135.047c-10.355,0-18.751-8.393-18.751-18.751c0-10.355,8.396-18.751,18.751-18.751h341.905
                              c10.355,0,18.751,8.396,18.751,18.751C495.702,260.641,487.307,269.033,476.951,269.033z"/>
                          </g>
                          <g>
                            <path d="M476.951,380.469H135.047c-10.355,0-18.751-8.393-18.751-18.751c0-10.355,8.396-18.751,18.751-18.751h341.905
                              c10.355,0,18.751,8.396,18.751,18.751C495.702,372.076,487.307,380.469,476.951,380.469z"/>
                          </g>
                          <g>
                            <path d="M278.154,491.906H135.047c-10.355,0-18.751-8.394-18.751-18.751c0-10.355,8.396-18.751,18.751-18.751h143.106
                              c10.355,0,18.751,8.396,18.751,18.751C296.905,483.512,288.509,491.906,278.154,491.906z"/>
                          </g>
                        </g>
                      </g>
                    </g>
                    </svg>
                  </button>
                </th>           


                <th>
                  <button class="btn bg-error border-error" @click="gotoProfile"><p class="px-2"> Delete </p>
                  <svg fill="#000" version="1.1" id="Capa_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" 
                    class="h-6 w-6" viewBox="0 0 468.36 468.36"
                    xml:space="preserve">
                  <g>
                    <g>
                      <path d="M381.048,64.229l-71.396,0.031L309.624,0L158.666,0.064l0.027,64.26l-71.405,0.031l0.024,60.056h293.76L381.048,64.229z
                        M189.274,30.652l89.759-0.04l0.016,33.66l-89.759,0.04L189.274,30.652z"/>
                      <path d="M87.312,468.36h293.76V139.71H87.312V468.36z M303.042,184.588h15.301v238.891h-15.301V184.588z M226.542,184.588h15.3
                        v238.891h-15.3V184.588z M150.042,184.588h15.3v238.891h-15.3V184.588z"/>
                    </g>
                  </g>
                  </svg>
                  </button>
                </th>
              </tr>

            </tbody>
          </table>
        </div>  
      </div>
      
    </div>  
  </div>
  <div id="igvSection" class="flex justify-center items-center h-144 w-full top-96">
    <div class="w-11/12 min-w-96 h-full bg-white shadow-md rounded-md">
      <div class="flex items-center justify-center w-full h-32 inherit">
        <h1 class="text-4xl font-bold link-hover"><a href="https://github.com/igvteam/igv.js/">IGV.js</a></h1>
      </div>
      <div class="flex items-center justify-center w-full h-96">
        <div class="skeleton w-5/6 h-full"></div>
      </div>
    </div>
    </div>

</div>

</template>

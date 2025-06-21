<script setup>
/* eslint-disable */
import { ScrollTo } from '@/scripts/OtherFunc';
import { backendLink } from '@/scripts/AccountFunc';
import { ref } from 'vue';
import { RouterLink } from 'vue-router';
import axios from "axios";

const uploadedFileRef = ref(33);

async function getStat()
{
  let response = await axios.get(backendLink+"/statistic")
      .catch(console.error);
  console.log(response)
  if (response) uploadedFileRef.value = response.data.files_uploaded;
  
}

getStat();

// TODO: Add a 404 page

</script>
<template>
  <div
    class="fixed -z-20 hero min-h-screen bg-base-200 bg-linear-to-br from-blue-100 to-purple-100"
  ></div>

  <div class="hero min-h-screen">
    <div class="flex justify-center items-center text-center align-items-top shadow-2xl py-4 md:py-8 lg:py-16 bg-white rounded-md max-w-11/12">
      <div class="max-w-6xl grid grid-cols-1">
        <div class="flex justify-center items-center w-full">
          <h1 class="text-xl sm:text-3xl md:text-4xl xl:text-5xl font-bold leading-normal w-11/12">
            HLA-SeqPredict: A Comprehensive ADR Risk Assessment Platform
          </h1>
        </div>
        <p class="py-6 px-8 text-sm sm:text-md md:text-lg">
            HLA-SeqPredict is a user-friendly web application designed to provide a comprehensive assessment of a patient's risk for developing adverse drug reactions (ADRs)<span class="hidden sm:inline">, particularly severe immune-mediated cutaneous ADRs like Stevens-Johnson Syndrome (SJS) and Toxic Epidermal Necrolysis (TEN).</span>
        </p>
        <div class="w-full h-16 flex justify-center items-center gap-10">
        <button @click="ScrollTo('infoPage')" class="btn btn-info shadow-sm btn-sm md:btn-md">More Info</button>
        <RouterLink to="signup" class="btn btn-accent shadow-sm btn-sm md:btn-md">Get Started</RouterLink>
        <a
            href="https://github.com/theme222/HLA-SeqPredict"
            class="btn btn-secondary btn-sm md:btn-md"
            >Github</a
        >
        </div>
      </div>
    </div>




  </div>




  <div class="w-full flex justify-center items-center min-h-screen" id="infoPage">
    <div class="w-11/12 max-w-300 bg-white rounded-md shadow-md">
      <h1 class="text-2xl md:text-xl lg:text-2xl leading-normal font-medium text-black text-center px-5 xl:px-20 md:py-10 py-8">
        This project <span class="hidden md:inline">introduces a web-based tool to support the prevention and management of HLA-related ADRs. It </span>contains 3 main components: 
      </h1>

      <div class="md:hidden grid grid-cols-1 gap-6 pb-5">

        <div class="flex justify-center items-center w-full ">
          <div class="w-11/12 h-full border-[1px] border-secondary rounded-md hover:shadow-xl transition-all duration-400">
            <div class="w-full flex justify-center">
              <div class="w-46/50 h-16 flex justify-between items-center">
                <RouterLink to="upload" class="text-2xl font-bold text-primary">HLA Genotyping Analysis</RouterLink>
                <img src="@/assets/analysis.svg" alt="" class="w-10 h-10">
              </div>
            </div>
            <p class="text-md px-5 pb-4">
              This component processes nucleotide sequence data from HLA-targeted sequencing and identifies HLA types, reporting their association with ADR risks.
            </p>
          </div>
        </div>

        <div class="flex justify-center items-center w-full ">
          <div class="w-11/12 h-full border-[1px] border-secondary rounded-md hover:shadow-xl transition-all duration-400">
            <div class="w-full flex justify-center">
              <div class="w-46/50 h-16 flex justify-between items-center">
                <RouterLink to="search" class="sm:text-2xl font-bold text-primary">Pharmacogenomic Knowledge Base</RouterLink>
                <img src="@/assets/sunburst.svg" alt="" class="w-10 h-10">
              </div>
            </div>
            <p class="text-sm sm:text-base px-5 pb-4">
              This component presents information on the relationships between each HLA type, associated drugs, and ADR risks in the form of a sunburst plot.This allows doctors and pharmacists to search for information and make informed decisions when prescribing medications.
            </p>
          </div>
        </div>
        
        <div class="flex justify-center items-center w-full ">
          <div class="w-11/12 h-full border-[1px] border-secondary rounded-md hover:shadow-xl transition-all duration-400">
            <div class="w-full flex justify-center">
              <div class="w-46/50 h-16 flex justify-between items-center">
                <RouterLink to="search" class="text-xl sm:text-2xl font-bold text-primary">HLA Prevalence Visualization</RouterLink>
                <img src="@/assets/histogram.svg" alt="" class="w-10 h-10">
              </div>
            </div>
            <p class="px-5 pb-4">
                This component presents data on the prevalence of each HLA type in the Thai population, including their association with ADR risks. This information supports policymakers, such as the NHSO, in considering expanding the coverage of HLA testing benefits.
            </p>
          </div>
        </div>
        
      </div>

      <div class="hidden md:flex h-144 justify-center items-center gap-6 xl:gap-12">
        <div class="w-15/50 group grid grid-cols-1 justify-center gap-10 max-h-133">
          <div class="border-secondary border-[1px] rounded-md group-hover:shadow-md h-95 min-h-max">
            <div class="h-20 w-full flex items-center">
              <img src="@/assets/analysis.svg" alt="" class="w-20 h-10">
            </div>
            <div class="w-full flex justify-center items-center">
              <RouterLink to="upload" class="w-11/12 font-bold text-xl lg:text-3xl text-primary cursor-pointer">HLA Genotyping Analysis</RouterLink>
            </div>          
            <div class="h-56 w-full flex justify-center items-center">
              <div class="w-11/12 lg:text-xl">
                This component processes nucleotide sequence data from HLA-targeted sequencing and identifies HLA types, reporting their association with ADR risks.
              </div>
            </div>
          </div>

          <div class="stats opacity-20 group-hover:border-2 border-primary group-hover:opacity-100 group-hover:shadow-xl transition-all duration-400 ease-in-out">
            <div class="stat">
              <div class="stat-title">Uploaded Files</div>
              <div class="stat-value">{{ uploadedFileRef }}</div>
              <div class="stat-desc">Successfully visualized with embedded IGV</div>
            </div>
          </div>
        </div>


        <div class="w-15/50 group grid grid-cols-1 justify-center gap-10">
          <div class="border-secondary border-[1px] rounded-md group-hover:shadow-md h-95 min-h-max">
            <div class="h-20 w-full flex items-center">
            <img src="@/assets/sunburst.svg" alt="sunburst chart" class="w-20 h-12">
          </div>
            <div class="w-full flex justify-center items-center">
              <RouterLink to="search" class="w-11/12 font-bold text-primary md:text-xl lg:text-2xl xl:text-3xl">Pharmacogenomic Knowledge Base</RouterLink>
            </div>          
            <div class="h-56 w-full flex justify-center items-center">
              <h2 class="w-11/12 text-sm md:text-base lg:text-xl xl:text-base"> 
                This component presents information on the relationships between each HLA type, associated drugs, and ADR risks in the form of a sunburst plot. <span class="hidden xl:inline">This allows doctors and pharmacists to search for information and make informed decisions when prescribing medications.</span>
              </h2>
            </div>
          </div>

        <div class="stats opacity-20 group-hover:border-2 border-primary  group-hover:opacity-100 group-hover:shadow-xl transition-all duration-400 ease-in-out ">
            <div class="stat">
              <div class="stat-title">Sunburst Entries</div>
              <div class="stat-value">175,363</div>
              <div class="stat-desc">Collected data from GeTH</div>
            </div>
        </div>
        </div>

        <div class="w-15/50 group grid grid-cols-1 justify-center gap-10">
          <div class="border-secondary border-[1px] rounded-md group-hover:shadow-md h-95 min-h-max">
            <div class="h-20 w-full flex items-center">
              <img src="@/assets/histogram.svg" alt="" class="w-20 h-12">
          </div>
            <div class="w-full flex justify-center items-center">
              <RouterLink to="search" class="w-11/12 font-bold text-primary text-xl lg:text-3xl">HLA Prevalence Visualization</RouterLink>
            </div>          
            <div class="h-56 w-full flex justify-center items-center">
              <h2 class="w-11/12 text-sm lg:text-base"> 
                This component presents data on the prevalence of each HLA type in the Thai population, including their association with ADR risks. This information supports policymakers, such as the NHSO, in considering expanding the coverage of HLA testing benefits.
              </h2>
            </div>
          </div>
        

        <div class="stats opacity-20 group-hover:border-2 border-primary  group-hover:opacity-100 group-hover:shadow-xl transition-all duration-400 ease-in-out ">
            <div class="stat">
              <div class="stat-title">Drug Listings</div>
              <div class="stat-value">18</div>
              <div class="stat-desc">Searchable using the search bar</div>
            </div>
        </div>
        </div>
      </div>
    </div>
  </div>  


</template>


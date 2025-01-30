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
</script>
<template>
  <div
    class="fixed -z-20 hero min-h-screen bg-base-200 bg-gradient-to-br from-blue-100 to-purple-100"
  ></div>

  <div class="hero min-h-screen">
    <div class="hero-content text-center align-items-top shadow-2xl py-16 bg-white rounded-md">
      <div class="max-w-6xl">
        <h1 class="text-5xl font-bold leading-normal">
          HLA-SeqPredict: A Comprehensive ADR Risk Assessment Platform
        </h1>
        <p class="py-6 px-8 text-lg">
            HLA-SeqPredict is a user-friendly web application designed to provide a comprehensive assessment of a patient's risk for developing adverse drug reactions (ADRs), particularly severe immune-mediated cutaneous ADRs like Stevens-Johnson Syndrome (SJS) and Toxic Epidermal Necrolysis (TEN).

        </p>
        <div class="w-full h-16  flex justify-center items-center gap-10">
        <button @click="ScrollTo('infoPage')" class="btn btn-info shadow">More Info</button>
        <RouterLink to="signup" class="btn btn-accent shadow">Get Started</RouterLink>
        <a
            href="https://github.com/theme222/HLA-SeqPredict"
            class="btn btn-secondary"
            >Github</a
        >
        </div>
      </div>
    </div>




  </div>




  <div class="w-full flex justify-center items-center min-h-screen">
    <div class="w-300 bg-white rounded-md" id="infoPage">
      <h1 class="text-2xl leading-normal font-medium text-black text-center h-32 px-20 py-10 flex justify-center items-center">
        This project introduces a web-based tool to support the prevention and management of HLA-related ADRs consisting of 3 main components: 
      </h1>

      <div class="flex h-144 justify-center items-center gap-12">

        <div class="w-14/50 group grid grid-cols-1 justify-center gap-10">
          <div class="border-secondary border-[1px] rounded-md group-hover:shadow-md">
            <div class="h-20 w-full flex items-center">
              <svg class="w-20 h-10" version="1.1" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" viewBox="0 0 121.7 122.88" style="enable-background:new 0 0 121.7 122.88" xml:space="preserve"><g><path d="M53.62,0c14.81,0,28.21,6,37.91,15.71c9.7,9.7,15.71,23.11,15.71,37.91c0,10.83-3.21,20.91-8.74,29.35l23.21,25.29 l-16,14.63L83.32,98.26c-8.5,5.67-18.72,8.98-29.7,8.98c-14.81,0-28.21-6-37.91-15.71C6,81.82,0,68.42,0,53.62 C0,38.81,6,25.41,15.7,15.7C25.41,6,38.81,0,53.62,0L53.62,0L53.62,0z M54.73,41.18l1.49-5.67l-8.46-2.97l17.18-9.21l9.14,17.42 l-8.99-2.8l-1.25,5.95L54.73,41.18L54.73,41.18z M83.18,81.99H29.67c-1.07,0-2.09-0.43-2.85-1.19c-0.72-0.72-1.17-1.73-1.17-2.83 L25.64,35.1h6.06v40.83h51.48V81.99L83.18,81.99z M69.92,53.32h7.75c0.44,0,0.79,0.36,0.79,0.79v16.14c0,0.44-0.36,0.79-0.79,0.79 h-7.75c-0.44,0-0.79-0.36-0.79-0.79V54.12C69.12,53.68,69.48,53.32,69.92,53.32L69.92,53.32L69.92,53.32z M39.38,43.44h7.75 c0.44,0,0.79,0.36,0.79,0.79v26.03c0,0.44-0.36,0.79-0.79,0.79h-7.75c-0.44,0-0.79-0.36-0.79-0.79V44.24 C38.59,43.8,38.94,43.44,39.38,43.44L39.38,43.44L39.38,43.44L39.38,43.44z M63.42,48.68v21.59c0,0.44-0.36,0.79-0.79,0.79h-7.74 c-0.44,0-0.79-0.36-0.79-0.79V46.01l8.38,2.49c0.13,0.04,0.26,0.07,0.39,0.1C63.04,48.63,63.23,48.66,63.42,48.68L63.42,48.68z M87.3,19.93C78.68,11.31,66.77,5.98,53.62,5.98c-13.15,0-25.06,5.33-33.68,13.95C11.31,28.55,5.98,40.46,5.98,53.62 c0,13.15,5.33,25.06,13.95,33.68c8.62,8.62,20.53,13.95,33.68,13.95c13.16,0,25.06-5.33,33.68-13.95 c8.62-8.62,13.95-20.53,13.95-33.68C101.25,40.46,95.92,28.55,87.3,19.93L87.3,19.93L87.3,19.93z"/></g></svg>
            </div>
            <div class="w-full flex justify-center items-center">
              <RouterLink to="upload" class="w-11/12 font-bold text-3xl text-primary cursor-pointer">HLA Genotyping Analysis</RouterLink>
            </div>          
            <div class="h-56 w-full flex justify-center items-center">
              <div class="w-11/12 text-xl">
                This component processes nucleotide sequence data from HLA-targeted sequencing and identifies HLA types, reporting their association with ADR risks.
              </div>
            </div>
          </div>

          <div class="stats opacity-20 group-hover:border-2 border-primary group-hover:opacity-100 group-hover:shadow-xl transition-all duration-400 ease-in-out ">
            <div class="stat">
              <div class="stat-title">Uploaded Files</div>
              <div class="stat-value">{{ uploadedFileRef }}</div>
              <div class="stat-desc">Successfully visualized with embedded IGV</div>
            </div>
          </div>
        </div>


        <div class="w-14/50 group grid grid-cols-1 justify-center gap-10">
          <div class="border-secondary border-[1px] rounded-md group-hover:shadow-md">
            <div class="h-20 w-full flex items-center">
              <svg class="w-20 h-12" fill="#000000" version="1.1" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" 
                viewBox="0 0 32 32" xml:space="preserve">
              <path id="chart--sunburst_1_" d="M16,31.36c-4.103,0-7.96-1.598-10.861-4.499c-0.141-0.141-0.141-0.369,0-0.51l3.289-3.288
                c-1.652-1.771-2.693-4.118-2.782-6.704H1c-0.199,0-0.36-0.161-0.36-0.36c0-4.103,1.597-7.959,4.498-10.86
                c0.027-0.027,0.059-0.05,0.093-0.067c0.185-0.093,0.227-0.107,3.704,3.357c1.77-1.653,4.119-2.694,6.705-2.783V1
                c0-0.199,0.161-0.36,0.36-0.36c8.47,0,15.36,6.89,15.36,15.36c0,0.199-0.161,0.36-0.36,0.36h-4.646
                c-0.187,5.427-4.566,9.807-9.994,9.994V31C16.36,31.199,16.199,31.36,16,31.36z M5.906,26.604c2.646,2.523,6.078,3.944,9.734,4.032
                v-4.282c-2.585-0.089-4.934-1.129-6.704-2.781L5.906,26.604z M16.36,21.348v4.286c5.031-0.186,9.088-4.243,9.274-9.274h-4.286
                C21.17,19.029,19.029,21.17,16.36,21.348z M9.446,23.063c1.64,1.522,3.807,2.482,6.194,2.57v-4.286
                c-1.208-0.08-2.307-0.562-3.166-1.313L9.446,23.063z M6.367,16.36c0.088,2.387,1.048,4.555,2.57,6.194l3.029-3.028
                c-0.751-0.858-1.233-1.958-1.314-3.166H6.367z M16,11.36c-2.559,0-4.64,2.082-4.64,4.64s2.082,4.64,4.64,4.64s4.64-2.081,4.64-4.64
                S18.559,11.36,16,11.36z M26.354,15.64h4.282C30.447,7.853,24.147,1.553,16.36,1.364v4.282C21.787,5.833,26.167,10.213,26.354,15.64
                z M21.348,15.64h4.286c-0.186-5.03-4.243-9.088-9.274-9.273v4.285C19.029,10.83,21.17,12.971,21.348,15.64z M6.367,15.64h4.285
                c0.081-1.208,0.563-2.307,1.314-3.166c-0.492-0.492-1.736-1.735-3.029-3.029C7.415,11.085,6.455,13.253,6.367,15.64z M1.364,15.64
                h4.282c0.089-2.585,1.13-4.934,2.782-6.703C7.15,7.66,5.915,6.424,5.396,5.906C2.873,8.553,1.452,11.984,1.364,15.64z M9.445,8.938
                c0.828,0.826,1.825,1.823,3.03,3.028c0.858-0.751,1.958-1.233,3.166-1.314V6.367C13.252,6.455,11.084,7.415,9.445,8.938z"/>
              <rect id="_Transparent_Rectangle" style="fill:none;" width="32" height="32"/>
              </svg>
          </div>
            <div class="w-full flex justify-center items-center">
              <RouterLink to="search" class="w-11/12 font-bold text-primary text-3xl">Pharmacogenomic Knowledge Base</RouterLink>
            </div>          
            <div class="h-56 w-full flex justify-center items-center">
              <h2 class="w-11/12 text-md"> 
                This component presents information on the relationships between each HLA type, associated drugs, and ADR risks in the form of a sunburst plot. This allows doctors and pharmacists to search for information and make informed decisions when prescribing medications.
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

        <div class="w-14/50 group grid grid-cols-1 justify-center gap-10">
          <div class="border-secondary border-[1px] rounded-md group-hover:shadow-md">
            <div class="h-20 w-full flex items-center">
              <svg class="w-20 h-12 cursor-pointer" xmlns="http://www.w3.org/2000/svg" viewBox="-0.5 -0.5 16 16" fill="none" stroke="#000000" stroke-linecap="round" stroke-linejoin="round" id="Chart-Histogram--Streamline-Tabler" height="16" width="16"><desc>Chart Histogram Streamline Icon: https://streamlinehq.com</desc><path d="M1.875 1.875v11.25h11.25" stroke-width="1"></path><path d="M12.5 11.25v1.875" stroke-width="1"></path><path d="M10 10v3.125" stroke-width="1"></path><path d="M7.5 8.125v5" stroke-width="1"></path><path d="M5 10v3.125" stroke-width="1"></path><path d="M1.875 6.875c3.75 0 3.125 -3.125 5.625 -3.125s1.875 3.125 5.625 3.125" stroke-width="1"></path></svg>
          </div>
            <div class="w-full flex justify-center items-center">
              <RouterLink to="search" class="w-11/12 font-bold text-primary text-3xl">HLA Prevalence Visualization</RouterLink>
            </div>          
            <div class="h-56 w-full flex justify-center items-center">
              <h2 class="w-11/12 text-md"> 
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


<script setup>
/* eslint-disable */
import { ref, onMounted} from "vue";
import { getAccountInfo } from "@/scripts/AccountFunc";
import { useRouter, useRoute } from "vue-router";

const router = useRouter();
const route = useRoute();

const error = ref();
const accountName = ref();
const accountEmail = ref();

const searchBarText = ref('');
const activeSuggestionList = ref([]);
let suggestionList;

const viewportSize = ref([0,0]);

async function SearchBarSetup()
{
  const histogramDataModule = await import("./data/histogram.json");
  const histogramData = histogramDataModule.default; // Access the JSON content

  const drugSet = new Set();
  for (let obj of Object.values(histogramData)) {
    for (let drug of obj.drug) {
      if (drug !== '') drugSet.add(drug);
    }
  }

  suggestionList = Array.from(drugSet);
}

onMounted(SearchBarSetup);

function GiveSuggestion()
{
  if (!suggestionList) return [];
  const compare = (word1, word2) => {
    // An implementation of the Wagner–Fischer algorithm

    if (word1.length > word2.length) word1 = word1.slice(0,word2.length);
    else if (word2.length > word1.length) word2 = word2.slice(0,word1.length);

    let n = word1.length + 1; // X axis
    let m = word2.length + 1; // Y axis
    const arr = Array.from({ length: n }, () => Array(m).fill(0));

    for (let i = 0; i < n; i++) arr[i][0] = i;
    for (let j = 0; j < m; j++) arr[0][j] = j;

    for (let i = 1; i < n; i++)
    {
      for (let j = 1; j < m; j++)
      {
        let subValue = Number(word1[i-1].toLowerCase() != word2[j-1].toLowerCase());
        arr[i][j] = Math.min(arr[i-1][j] + 1, arr[i][j-1] + 1, arr[i-1][j-1] + subValue);
      }
    }
    return arr[n-1][m-1];
  }
  const output = [];
  for (let index in suggestionList)
  {
    output.push([compare(searchBarText.value,suggestionList[index]), suggestionList[index]]);
  }
  activeSuggestionList.value = output.sort((a,b) => a[0] - b[0]).slice(0,3);
}

function DoSearchQuery(_query) {
  if (suggestionList && suggestionList.includes(_query)) searchBarText.value = _query
  else if (_query != '')  searchBarText.value = activeSuggestionList.value[0][1];
  router.push({name: "search", query: {q: searchBarText.value || searchBarText.value}})
  const searchEvent = new CustomEvent("searchEvent", {detail: searchBarText.value || searchBarText.value});
  setTimeout(() => window.dispatchEvent(searchEvent),500)
}

getAccountInfo().then((data) => {
  if (data) {
    accountName.value = data["name"];
    accountEmail.value = data["email"];
  }
});

function ClearError() {
  window.sharedData.latestError = null;
  error.value = null;
}

setInterval(() => {
  if (!window.sharedData.latestError) 
  {
    error.value = null;
    return;
  }
  else
  {
    if (window.sharedData.latestError.response == null) 
    {
      error.value = window.sharedData.latestError.message;
      return;
    }
    error.value = `${window.sharedData.latestError.response.data} : ${window.sharedData.latestError.response.status}`;
  }
}, 200);

window.addEventListener("resize", () =>
{
  viewportSize.value = [window.innerWidth, window.innerHeight];
})



</script>

<template>
<!-- Navbar -->
<div class="z-69420 navbar flex justify-between bg-primary fixed">
    <div class="flex sm:min-w-54 xl:min-w-88 justify-start items-center">
      <RouterLink class="btn btn-ghost text-2xl text-white font-bold" to="home">HLA-SeqPredict</RouterLink>
    </div>
    <div class="justify-center items-center w-4/5 max-w-200 hidden sm:flex">
      <div class="w-full h-12 flex gap-1 px-5 items-center">
        <label class="group flex input w-full items-center justify-between px-5 dropdown hover:dropdown-open dropdown-bottom dropdown-end">
          <div class="w-full flex items-center">
            <input
              type="text"
              id="searchBar"
              name="searchBar"
              class="text-black px-2 w-full"
              placeholder="Search drug name"
              v-model="searchBarText"
              @keyup="GiveSuggestion"
              @keyup.enter="DoSearchQuery(searchBarText)"
            />
          </div>
          <ul class="dropdown-content w-full grid grid-cols-1 justify-center gap-0.5">
            <li class="h-4"></li>
            <li v-for="i in activeSuggestionList" :key='i[1]' @click="DoSearchQuery(i[1])" class="shadow-lg h-12 w-full flex justify-normal items-center btn bg-white hover:bg-slate-900 hover:text-white">
              <p class="text-lg px-2"> {{ i[1] }} </p>
            </li>
          </ul>
        </label>
        <div class="btn btn-neutral flex justify-center items-center w-14" @click="DoSearchQuery(searchBarText)">
          <img src="@/assets/search.svg" alt="search" class="w-6 h-6">
        </div>
      </div>
    </div>
    <div class="flex min-w-35 lg:min-w-88 justify-end">
      <ul class="menu menu-horizontal w-full px-1">
        <li class="hidden lg:block">
          <RouterLink to="dashboard" class="font-bold text-white">Dashboard</RouterLink>
        </li>
        <li class="hidden lg:block">
          <RouterLink to="upload" class="font-bold text-white">Upload Files</RouterLink>
        </li>
        <li>
          <details>
            <summary class="text-white font-bold">
              {{ accountName || "Account" }}
            </summary>
            <ul class="p-2 shadow-sm menu dropdown-content bg-base-100 rounded-box w-28 right-0 bg-ghost">
              <li class="lg:hidden"><RouterLink to="dashboard">Dashboard</RouterLink></li>
              <li class="lg:hidden"><RouterLink to="upload">Upload</RouterLink></li>
              <li class="sm:hidden"><RouterLink to="search">Search</RouterLink></li>
              <li><RouterLink to="profile">Profile</RouterLink></li>
              <li><RouterLink to="signup">Sign up</RouterLink></li>
              <li><RouterLink to="login">Login</RouterLink></li>
            </ul>
          </details>
        </li>
      </ul>
    </div>
</div>
<!-- Navbar -->

<!-- Background -->
<div
  class="fixed -z-20 hero min-h-screen bg-base-200 bg-linear-to-br from-blue-50 to-purple-50"
></div>
<!-- Background -->

<router-view />

<!-- Error View -->
<div
  class="bottom-5 w-full flex justify-end items-center fixed z-69420"
  v-if="error"
>
  <div class="px-5">
    <div role="alert" class="alert alert-error">
      <svg
        xmlns="http://www.w3.org/2000/svg"
        class="h-6 w-6 shrink-0 stroke-current"
        @click="ClearError"
        fill="none"
        viewBox="0 0 24 24"
      >
        <path
          stroke-linecap="round"
          stroke-linejoin="round"
          stroke-width="2"
          d="M10 14l2-2m0 0l2-2m-2 2l-2-2m2 2l2 2m7-2a9 9 0 11-18 0 9 9 0 0118 0z"
        />
      </svg>
      <span>{{ error }}</span>
    </div>
  </div>
</div>
<!-- Error View -->

<!-- Screen size debug -->
<div class="p-4 text-center text-white bg-info rounded-lg text-lg font-bold fixed bottom-0 z-9001" hidden>
    <span class="block sm:hidden">Extra Small (xs)</span>
    <span class="hidden sm:block md:hidden">Small (sm)</span>
    <span class="hidden md:block lg:hidden">Medium (md)</span>
    <span class="hidden lg:block xl:hidden">Large (lg)</span>
    <span class="hidden xl:block">Extra Large (xl)</span>
    <span>{{ viewportSize[0] }}x{{ viewportSize[1] }}</span>
</div>
<!-- Screen size debug -->
</template>


<script setup>
/* eslint-disable */
import { ref } from "vue";
import { getAccountInfo } from "@/scripts/AccountFunc";
import { useRouter, useRoute } from 'vue-router'

const router = useRouter()
const route = useRoute()

const error = ref();
const accountName = ref();
const accountEmail = ref();

getAccountInfo().then((data) => {
  if (data)
  {
  accountName.value = data['name'];
  accountEmail.value = data['email']
  }
});

function ClearError()
{
  window.sharedData.latestError = null
}

setInterval(() => {error.value = window.sharedData.latestError}, 100)

</script>

<template>
  
<div class="z-10 navbar bg-base-100 bg-primary fixed">
  <div class="flex-1">
    <RouterLink class="btn btn-ghost text-2xl text-white font-bold" to="home">ADR-Prediction</RouterLink>
  </div>
  <div class="flex-none">
    <ul class="menu menu-horizontal px-1">
      <li><RouterLink to="dashboard" class="font-bold text-white">Dashboard</RouterLink></li>
      <li><RouterLink to="upload" class="font-bold text-white">Upload Files</RouterLink></li>
      <li>
        <details>
          <summary class="text-white font-bold">
            {{ accountName || "Account" }}
          </summary>
          <ul class="p-2 shadow menu dropdown-content bg-base-100 rounded-box w-28 right-0 bg-ghost">
            <li><RouterLink to="profile">Profile</RouterLink></li>
            <li><RouterLink to="signup">Sign up</RouterLink></li>
            <li><RouterLink to="login">Login</RouterLink></li>
          </ul>
        </details>
      </li>
    </ul>
  </div>
</div>
  <div class="fixed -z-20 hero min-h-screen bg-base-200 bg-gradient-to-br  from-blue-50 to-purple-50"></div>

  <router-view />

  <div class="bottom-5 w-full flex justify-end items-center fixed z-50" v-if="error">
    <div class="px-5">
      <div role="alert" class="alert alert-error">
          <svg
            xmlns="http://www.w3.org/2000/svg"
            class="h-6 w-6 shrink-0 stroke-current"
            @click="ClearError"
            fill="none"
            viewBox="0 0 24 24">
            <path
              stroke-linecap="round"
              stroke-linejoin="round"
              stroke-width="2"
              d="M10 14l2-2m0 0l2-2m-2 2l-2-2m2 2l2 2m7-2a9 9 0 11-18 0 9 9 0 0118 0z" />
          </svg>
          <span>{{ error }}</span>
        </div>
    </div>
  </div>
  
</template>


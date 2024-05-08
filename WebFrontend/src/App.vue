<script setup>
/* eslint-disable */
import { ref } from "vue";
import { getAccountInfo } from "@/scripts/AccountFunc";
import { useRouter, useRoute } from 'vue-router'

const router = useRouter()
const route = useRoute()

const accountName = ref();
const accountEmail = ref();

getAccountInfo().then((data) => {
  accountName.value = data['name'];
  accountEmail.value = data['email']
});

function redirectToProfile(){
  if (accountName.value){
    router.push({name:"profile"})
  }
}

</script>

<template>
  <div class="background"></div>

  <div class="titleWrapper">
    <h1 class="titleText">
      <RouterLink to="home">HLA-ADR Prediction</RouterLink>
    </h1>

    <div to="profile" class="accountName" @click="redirectToProfile">{{ accountName || "No account" }}</div>
  </div>
  <router-view />
</template>

<style>
.background{
  position: fixed;
  top: 0px;
  left: 0px;
  padding: 100%;
  background-color: #EBEBEB;
}

.titleWrapper{
    font-family: Arial, Helvetica, sans-serif;
    position: fixed;
    top: 0px;
    left: 0px;
    background: #FF6700 ;
    width: 100%;
    padding-top: 1%;
    padding-bottom: 1%;
    z-index: 69420;
}

.titleText{
    display: inline;
    padding-left: 1%;
    color: white;
    font-weight: 800;
}



</style>
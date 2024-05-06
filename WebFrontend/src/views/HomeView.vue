<script setup>
/* eslint-disable */
let func = "This is so that autocomplete doesn't do some bullshit";
import { ref } from "vue";
import { reactive } from "vue";
import axios from "axios";
import Cookies from "js-cookie";
import { useRouter, useRoute } from "vue-router";
import { getAccountInfo } from "@/scripts/AccountFunc";

const router = useRouter();
const route = useRoute();

// Hash function
import { sha256 } from "js-sha256";

function sendToLogin() {
  router.push({ name: "login" });
}

function sendToSignup() {
  router.push({ name: "signup" });
}

const accountName = ref(null);
getAccountInfo().then((name) => {
  accountName.value = name;
});
</script>

<template>
  <div class="background"></div>

  <div class="titleWrapper">
    <h1 class="titleText">
      <RouterLink to="home">HLA-ADR Prediction</RouterLink>
    </h1>

    <div class="accountName">{{ accountName || "No account" }}</div>
  </div>

  <h1 class="heading">
    A prediction of adverse drug reactions through processing genotype
  </h1>

  <div class="description">
    This project is created by Sira Tongsima to help provide a quick and easy
    way to process raw sequencing data check out the source code at my
    <a href="https://github.com/theme222/ADR-Prediction" style="background: #36382e; padding: 0.3%; border-radius: 5px;">github</a>
  </div>

  <div class="info">Please start by</div>

  <RouterLink to="login" class="redirectButton" style="left: 30%"
    >Login</RouterLink
  >
  <div
    style="
      position: absolute;
      top: 32%;
      font-size: 1.7em;
      font-weight: bold;
      left: 50%;
      transform: translate(-50%);
    "
  >
    or
  </div>
  <RouterLink to="signup" class="redirectButton" style="left: 70.5%"
    >Signup</RouterLink
  >
  <div
    style="
      position: absolute;
      top: 39%;
      font-size: 1.1em;
      font-weight: bold;
      left: 50%;
      transform: translate(-50%);
    "
  >
    have an account?
  </div>
  <RouterLink to="upload" class="redirectButton" style="left: 30%; top: 42%"
    >Upload</RouterLink
  >
  <RouterLink to="dashboard" class="redirectButton" style="left: 70.6%; top: 42%"
    >Dashboard</RouterLink
  >
</template>

<style>
.info {
  position: absolute;
  top: 25%;
  left: 50%;
  transform: translate(-50%);
  font-weight: 900;
  font-size: 1.5em;
}

.background {
  position: fixed;
  top: 0px;
  left: 0px;
  padding: 100%;
  background-color: #ebebeb;
}

.titleWrapper {
  z-index: 983458934509834598;
  position: fixed;
  top: 0px;
  left: 0px;
  height: 40px;
  background: #ff6700;
  width: 100%;
  padding-top: 1.5%;
  padding-bottom: 1%;
}

.titleText {
  font-family: Arial, Helvetica, sans-serif;
  position: absolute;
  display: inline;
  color: white;
  font-weight: 800;
  transform: translate(0, -50%);
}

.loginLink {
  position: inherit;
  top: 2.5%;
  left: 85%;
}

.signupLink {
  position: inherit;
  top: 2.5%;
  left: 92%;
}

.heading {
  font-family: Arial, Helvetica, sans-serif;
  position: absolute;
  white-space: nowrap;
  top: 8%;
  left: 50%;
  transform: translate(-50%);
  overflow: unset;
  font-weight: 900;
  font-size: 2.3em;
}

/* The submit button */
.redirectButton {
  position: absolute;
  top: 30%;
  left: 50%;
  font-size: 1.2em;
  display: inline-block;
  outline: none;
  cursor: pointer;
  font-weight: 600;
  border-radius: 3px;
  padding: 40px 80px;
  border: 0;
  color: #3a4149;
  background: #ff6700;
  line-height: 1.15;
  transform: translate(-50%);
}

.redirectButton:hover {
  transition: all 0.1s ease;
  box-shadow: 0 0 0 0 #fff, 0 0 0 3px #1de9b6;
}

.description {
  position: absolute;
  left: 50%;
  top: 18%;
  text-align: center;
  transform: translate(-50%);
}
</style>
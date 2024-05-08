<script setup>
/* eslint-disable */
let func = "This is so that autocomplete doesn't do some bullshit";
import { ref } from "vue";
import { reactive } from "vue";
import axios from "axios";
import { useRouter, useRoute } from "vue-router";
import { backendLink } from "@/scripts/AccountFunc";

const router = useRouter();
const route = useRoute();

// Hash function
import { sha256 } from "js-sha256";
import Cookies from "js-cookie";

const formData = ref({
  name: "",
  email: "",
  password: "",
  confirmPassword: "",
});

const formCheck = ref({
  // false - doesn't exist, true - Good
  name: true,
  // 0 - no @, 1 - Good, 2 - Already Exists
  email: true,
  // false - short, true - good
  password: true,
  // false - not the same, true - good
  confirmPassword: true,
});

function checkValidName() {
  formCheck.value.name = formData.value.name.length > 0;
}

function checkValidEmail() {
  if (!formData.value.email.includes("@")) {
    formCheck.value.email = 0;
  } else {
    formCheck.value.email = 1;
    axios
      .post(backendLink+"/api/checkEmail", {
        email: formData.value.email,
      })
      .then((response) => {
        if (response.data["value"]) {
          formCheck.value.email = 2;
        }
      })
      .catch((error) => {
        alert("Server offline");
        console.log(error);
      });
  }
}

function checkValidPassword() {
  formCheck.value.password = formData.value.password.length >= 6;
  formCheck.value.confirmPassword =
    formData.value.password == formData.value.confirmPassword;
}

function allOfTheAbove() {
  checkValidName();
  checkValidEmail();
  checkValidPassword();

  //exampleData()

  let a = 0;
  for (let key in formCheck.value) {
    if (formCheck.value[key]) {
      a += 1;
    }
  }
  if (a == 4) {
    axios
      .post(backendLink+"/api/signup", {
        name: formData.value.name,
        email: formData.value.email,
        password: sha256(formData.value.password),
      })
      .then((response) => {
        router.push({ name: "login" });
      })
      .catch((error) => {
        console.error("Error signing up user:", error);
        // Handle error (e.g., display error message to user)
      });
  }
}
</script>

<template>


  <h2 class="signupHeading">Create an account</h2>

  <div class="inputForm" action="">
    <div>
      <input
        maxlength="100"
        type="text"
        id="name"
        class="box"
        v-model="formData.name"
        v-if="formCheck.name"
        required
        placeholder="Name / Username"
        @change="checkValidName"
      />
      <div v-else style="color: red; font-size: 0.9em">
        <input
          maxlength="100"
          type="text"
          id="name"
          class="youDidSomethingBadBox"
          v-model="formData.name"
          required
          placeholder="Name / Username"
          @change="checkValidName"
        />
        Please input a name / username
      </div>
    </div>
    <br />
    <div>
      <input
        maxlength="100"
        type="email"
        id="email"
        class="box"
        v-model="formData.email"
        v-if="formCheck.email == 1"
        required
        placeholder="Email"
        @change="checkValidEmail"
      />
      <div v-else style="color: red; font-size: 0.9em">
        <input
          maxlength="100"
          type="email"
          id="email"
          class="youDidSomethingBadBox"
          v-model="formData.email"
          required
          placeholder="Email"
          @change="checkValidEmail"
        />
        <div v-if="formCheck.email == 0">Please input a valid email</div>
        <div v-else>You can only register one account per email</div>
      </div>
    </div>
    <br />
    <div>
      <input
        maxlength="100"
        type="password"
        id="password"
        class="box"
        v-model="formData.password"
        v-if="formCheck.password"
        required
        placeholder="Password"
        @change="checkValidPassword"
      />
      <div v-else style="color: red; font-size: 0.9em">
        <input
          maxlength="100"
          type="password"
          id="password"
          class="youDidSomethingBadBox"
          v-model="formData.password"
          required
          placeholder="Password"
          @change="checkValidPassword"
        />
        Password should be atleast 6 characters
      </div>
    </div>
    <br />
    <div>
      <input
        maxlength="100"
        type="password"
        id="confirmPassword"
        class="box"
        v-model="formData.confirmPassword"
        v-if="formCheck.confirmPassword"
        required
        placeholder="Confirm password"
        @change="checkValidPassword"
      />
      <div v-else style="color: red; font-size: 0.9em">
        <input
          maxlength="100"
          type="password"
          id="password"
          class="youDidSomethingBadBox"
          v-model="formData.confirmPassword"
          required
          placeholder="Confirm password"
          @change="checkValidPassword"
        />
        Password does not match
      </div>
    </div>

    <button @click="allOfTheAbove" class="signupSubmitButton">Sign Up</button>
    <RouterLink to="login" class="redirectToLogin">
      Have an account? <u>Login</u></RouterLink
    >
  </div>
</template>

<style>
a:link {
  color: white;
  text-decoration: none;
}

a:visited {
  color: white;
  text-decoration: none;
}

a:hover {
  color: white;
}

a:active {
  color: white;
  text-decoration: none;
}

.background {
  position: fixed;
  top: 0px;
  left: 0px;
  padding: 100%;
  background-color: #ebebeb;
}

.signupHeading {
  font-family: Arial, Helvetica, sans-serif;
  position: absolute;
  top: 6%;
  left: 50%;
  transform: translate(-50%);
  font-weight: 900;
  font-size: 2.6em;
}

/* The input box thing lol */

.box {
  border: 3px solid #000;
  border-radius: 5px;
  height: 60px;
  line-height: normal;
  color: #282828;
  display: block;
  width: 500px;
  box-sizing: border-box;
  user-select: auto;
  font-size: 16px;
  padding: 0 6px;
  padding-left: 12px;
}
.box:focus {
  border: 3px solid #2b4570;
}

.youDidSomethingBadBox {
  border: 3px solid red;
  border-radius: 5px;
  height: 60px;
  line-height: normal;
  color: #282828;
  display: block;
  width: 500px;
  box-sizing: border-box;
  user-select: auto;
  font-size: 16px;
  padding: 0 6px;
  padding-left: 12px;
}

/* The submit button */
.signupSubmitButton {
  position: inherit;
  top: 130%;
  left: 50%;
  display: inline-block;
  outline: none;
  cursor: pointer;
  font-weight: 600;
  border-radius: 3px;
  padding: 12px 24px;
  border: 0;
  color: #3a4149;
  background: #ff6700;
  line-height: 1.15;
  font-size: 1em;
  transform: translate(-50%);
}
.signupSubmitButton:hover {
  transition: all 0.1s ease;
  box-shadow: 0 0 0 0 #fff, 0 0 0 3px #1de9b6;
}

.redirectToLogin {
  color: white;
  position: absolute;
  top: 110%;
  left: 50%;
  transform: translate(-50%);
  background: #ff6700;
  padding: 5px;
  border-radius: 3px;
}
</style>
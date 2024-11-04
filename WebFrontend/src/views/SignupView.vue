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
    formCheck.value.email = false;
  } else {
    formCheck.value.email = true;
    axios
      .post(backendLink + "/api/checkEmail", {
        email: formData.value.email,
      })
      .then((response) => {
        if (response.data["value"]) {
          formCheck.value.email = false;
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
      .post(backendLink + "/api/signup", {
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

  <div class="absolute w-full h-144 top-32 flex justify-center items-center">
    <div class="w-6/12 min-w-96 h-full shadow-md bg-white rounded-md">
      
      <div class="flex justify-center items-center w-full h-16 relative top-8 ">
        <h1 class="font-bold text-4xl">Create an account</h1>
      </div>
         
      <div class="relative flex justify-center items-center w-full h-80 top-16">
        <div class="grid gap-0 grid-rows-4 h-full w-10/12"> 

          
        <div v-if="formCheck.name">
          <label class="input input-bordered input-primary flex items-center gap-2">
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" class="w-4 h-4 opacity-70"><path d="M8 8a3 3 0 1 0 0-6 3 3 0 0 0 0 6ZM12.735 14c.618 0 1.093-.561.872-1.139a6.002 6.002 0 0 0-11.215 0c-.22.578.254 1.139.872 1.139h9.47Z" /></svg>          
            <input v-model="formData.name" type="text" class="grow" placeholder="Name" @change="checkValidName" />
          </label>
        </div>
        <div v-else>
          <label class="input input-bordered input-error flex items-center gap-2">
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" class="w-4 h-4 opacity-70"><path d="M8 8a3 3 0 1 0 0-6 3 3 0 0 0 0 6ZM12.735 14c.618 0 1.093-.561.872-1.139a6.002 6.002 0 0 0-11.215 0c-.22.578.254 1.139.872 1.139h9.47Z" /></svg>          
            <input v-model="formData.name" type="text" class="grow" placeholder="Name" @change="checkValidName" />
          </label>
          <div class="label"> <span class="label-text-alt text-error">Please input a name / username</span> </div>
        </div>
        
        <div v-if="formCheck.email">
          <label class="input input-bordered input-primary flex items-center gap-2">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" class="w-4 h-4 opacity-70"><path d="M2.5 3A1.5 1.5 0 0 0 1 4.5v.793c.026.009.051.02.076.032L7.674 8.51c.206.1.446.1.652 0l6.598-3.185A.755.755 0 0 1 15 5.293V4.5A1.5 1.5 0 0 0 13.5 3h-11Z" /><path d="M15 6.954 8.978 9.86a2.25 2.25 0 0 1-1.956 0L1 6.954V11.5A1.5 1.5 0 0 0 2.5 13h11a1.5 1.5 0 0 0 1.5-1.5V6.954Z" /></svg>            <input v-model="formData.email" type="text" class="grow" placeholder="Email" @change="checkValidEmail" />
          </label>
        </div>
        <div v-else>
          <label class="input input-bordered input-error flex items-center gap-2">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" class="w-4 h-4 opacity-70"><path d="M2.5 3A1.5 1.5 0 0 0 1 4.5v.793c.026.009.051.02.076.032L7.674 8.51c.206.1.446.1.652 0l6.598-3.185A.755.755 0 0 1 15 5.293V4.5A1.5 1.5 0 0 0 13.5 3h-11Z" /><path d="M15 6.954 8.978 9.86a2.25 2.25 0 0 1-1.956 0L1 6.954V11.5A1.5 1.5 0 0 0 2.5 13h11a1.5 1.5 0 0 0 1.5-1.5V6.954Z" /></svg>            <input v-model="formData.email" type="text" class="grow" placeholder="Email" @change="checkValidEmail" />
          </label>
          <div class="label"> <span class="label-text-alt text-error">Email is either invalid or already being used</span> </div>
        </div>

        <div v-if="formCheck.password">
          <label class="input input-bordered input-primary flex items-center gap-2">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" class="w-4 h-4 opacity-70"><path fill-rule="evenodd" d="M14 6a4 4 0 0 1-4.899 3.899l-1.955 1.955a.5.5 0 0 1-.353.146H5v1.5a.5.5 0 0 1-.5.5h-2a.5.5 0 0 1-.5-.5v-2.293a.5.5 0 0 1 .146-.353l3.955-3.955A4 4 0 1 1 14 6Zm-4-2a.75.75 0 0 0 0 1.5.5.5 0 0 1 .5.5.75.75 0 0 0 1.5 0 2 2 0 0 0-2-2Z" clip-rule="evenodd" /></svg>          
          <input v-model="formData.password" type="password" class="grow" placeholder="Password" @change="checkValidPassword" />
          </label>
        </div>
        <div v-else>
          <label class="input input-bordered input-error flex items-center gap-2">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" class="w-4 h-4 opacity-70"><path fill-rule="evenodd" d="M14 6a4 4 0 0 1-4.899 3.899l-1.955 1.955a.5.5 0 0 1-.353.146H5v1.5a.5.5 0 0 1-.5.5h-2a.5.5 0 0 1-.5-.5v-2.293a.5.5 0 0 1 .146-.353l3.955-3.955A4 4 0 1 1 14 6Zm-4-2a.75.75 0 0 0 0 1.5.5.5 0 0 1 .5.5.75.75 0 0 0 1.5 0 2 2 0 0 0-2-2Z" clip-rule="evenodd" /></svg>          
          <input v-model="formData.password" type="password" class="grow" placeholder="Password" @change="checkValidPassword" />
          </label>
          <div class="label"> <span class="label-text-alt text-error">Password is must be at least 6 characters long</span> </div>
        </div>

        <div v-if="formCheck.confirmPassword">
          <label class="input input-bordered input-primary flex items-center gap-2">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" class="w-4 h-4 opacity-70"><path fill-rule="evenodd" d="M14 6a4 4 0 0 1-4.899 3.899l-1.955 1.955a.5.5 0 0 1-.353.146H5v1.5a.5.5 0 0 1-.5.5h-2a.5.5 0 0 1-.5-.5v-2.293a.5.5 0 0 1 .146-.353l3.955-3.955A4 4 0 1 1 14 6Zm-4-2a.75.75 0 0 0 0 1.5.5.5 0 0 1 .5.5.75.75 0 0 0 1.5 0 2 2 0 0 0-2-2Z" clip-rule="evenodd" /></svg>          
          <input v-model="formData.confirmPassword" type="password" class="grow" placeholder="Confirm Password" @change="checkValidPassword" />
          </label>
        </div>
        <div v-else>
          <label class="input input-bordered input-error flex items-center gap-2">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="currentColor" class="w-4 h-4 opacity-70"><path fill-rule="evenodd" d="M14 6a4 4 0 0 1-4.899 3.899l-1.955 1.955a.5.5 0 0 1-.353.146H5v1.5a.5.5 0 0 1-.5.5h-2a.5.5 0 0 1-.5-.5v-2.293a.5.5 0 0 1 .146-.353l3.955-3.955A4 4 0 1 1 14 6Zm-4-2a.75.75 0 0 0 0 1.5.5.5 0 0 1 .5.5.75.75 0 0 0 1.5 0 2 2 0 0 0-2-2Z" clip-rule="evenodd" /></svg>          
          <input v-model="formData.confirmPassword" type="password" class="grow" placeholder="Confirm Password" @change="checkValidPassword" />
          </label>
          <div class="label"> <span class="label-text-alt text-error">Password does not match</span> </div>
        </div>

        </div>
      </div>

      <div class="relative grid grid-cols-1 place-content-center w-full h-8 top-12 text-secondary text-xs"> 
        <p class="w-full text-center py-2">
        Have an account? <RouterLink to="login" class="link">Login</RouterLink>.
        </p>
        <p class="w-full text-center">
        By clicking “Sign up”, you agree to the <RouterLink to="pdpa" class="link">Personal Data Protection Act</RouterLink>.
        </p>
      </div>

      <div class="relative flex justify-center items-center w-full h-16 top-20"> 
        <button @click="allOfTheAbove" type="submit" class="btn min-w-32 h-5/6 btn-primary text-white">Sign up</button>
      </div>



    </div>
  </div>

  
</template>


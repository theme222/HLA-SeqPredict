<script setup>
/* eslint-disable */
import { ref } from 'vue'
import axios from 'axios'
import Cookies from 'js-cookie'
import { useRouter, useRoute } from 'vue-router'
import { backendLink } from '@/scripts/AccountFunc'

const router = useRouter()
const route = useRoute()

// Hash function
import { sha256 } from 'js-sha256'

const formData = ref({
  email: '',
  password: ''
})


const formCheck = ref(true)

function GetExpirationDate(duration){
  const expirationDate = new Date()
  expirationDate.setTime(expirationDate.getTime() + duration*1000)
  return expirationDate
}

function sendLoginData(){
  axios.post(backendLink+'/api/login', {email : formData.value.email, password : sha256(formData.value.password)})
    .then(response => {
      if (response.data["success"]){
        formCheck.value = true
        Cookies.set('session_cookie', response.data["sessionCookie"], {expires : GetExpirationDate(response.data['duration'])})
        router.push({name: 'upload'})
        // REDIRECT TO MAIN PAGE :)

      } else {
        formCheck.value = false
      }
    })
    .catch(error => {
      console.error('Error signing up user:', error);
      // Handle error (e.g., display error message to user)
    });
}

</script>

<template>

    <div class="absolute w-full h-96 top-32 flex justify-center items-center">
    <div class="w-6/12 min-w-96 h-full shadow-md bg-white rounded-md">
      
      <div class="flex justify-center items-center w-full h-16 relative top-4 ">
        <h1 class="font-bold text-4xl">Login</h1>
      </div>
         
      <div class="relative flex justify-center items-center w-full h-32 top-10">
        <div class="grid gap-20 grid-rows-4 h-full w-10/12"> 
        
        <div v-if="formCheck">
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

        <div v-if="formCheck">
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

        </div>
      </div>

      <div class="relative grid grid-cols-1 place-content-center w-full h-8 top-12 text-secondary text-xs"> 
        <p class="w-full text-center">
        Don't have an account? <RouterLink to="signup" class="link">Sign up</RouterLink>.
        </p>
      </div>

      <div class="relative flex justify-center items-center w-full h-16 top-20"> 
        <button @click="sendLoginData" type="submit" class="btn min-w-32 h-5/6 btn-primary text-white">Login</button>
      </div>



    </div>
  </div>
    
</template>


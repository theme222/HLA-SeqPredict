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

  
    <h2 class="loginHeading">
      Login
    </h2>

    <div class="inputForm" action="">
      <div>      

        <input maxlength="100" type="email" id="email" class="box" v-model="formData.email" v-if="formCheck" required placeholder="Email">
        <div v-else style="color: red; font-size: 0.9em;">
          <input maxlength="100" type="email" id="email" class="youDidSomethingBadBox" v-model="formData.email" required placeholder="Email"> 
        </div>
        
      </div>
      <br>
      <div>

        <input maxlength="100" type="password" id="password" class="box" v-model="formData.password" v-if="formCheck" required placeholder="Password">
        <div v-else style="color: red; font-size: 0.9em;">
          <input maxlength="100" type="password" id="password" class="youDidSomethingBadBox" v-model="formData.password" required placeholder="Password"> 
          Email or password is incorrect
        </div>

      </div>

      <button @click="sendLoginData" class="loginSubmitButton">Login</button>
      <RouterLink to="signup" class="redirectToSignup" > Don't have an account? <u>Sign up</u></RouterLink>
      
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

.background{
  position: fixed;
  top: 0px;
  left: 0px;
  padding: 100%;
  background-color: #EBEBEB;
}


.loginLink{
  position: inherit;
  top: 2.5%;
  left: 85%;
}

.signupLink{
  position: inherit;
  top: 2.5%;
  left: 92%;
}

.loginHeading{
  font-family: Arial, Helvetica, sans-serif;
  position: absolute;
  top: 6%;
  left: 50%;
  transform: translate(-50%);
  font-weight: 900;
  font-size: 2.6em;
}

.inputForm{
  position: absolute;
  top: 15%;
  left: 50%;
  transform: translate(-50%);
}

/* The input box thing lol */ 

.box{
  border: 3px solid #000;
  border-radius: 5px;
  height: 60px;
  line-height: normal;
  color: #282828;
  display: block;
  width: 500px;
  box-sizing:border-box;
  user-select: auto;
  font-size: 16px;
  padding: 0 6px;
  padding-left: 12px;
}
.box:focus{
  border: 3px solid #2B4570;
}

.youDidSomethingBadBox {
  border: 3px solid red;
  border-radius: 5px;
  height: 60px;
  line-height: normal;
  color: #282828;
  display: block;
  width: 500px;
  box-sizing:border-box;
  user-select: auto;
  font-size: 16px;
  padding: 0 6px;
  padding-left: 12px;
}                
                
/* The submit button */
.loginSubmitButton{
  position: inherit;
  top: 118%;
  left: 50%;
  display: inline-block;
  outline: none;
  cursor: pointer;
  font-weight: 600;
  border-radius: 3px;
  padding: 12px 24px;
  border: 0;
  color: #3a4149;
  background: #FF6700;
  line-height: 1.15;
  font-size: 1em;
  transform: translate(-50%);
}
  .loginSubmitButton:hover {
  transition: all .1s ease;
  box-shadow: 0 0 0 0 #fff, 0 0 0 3px #1de9b6;
}

.redirectToSignup {
  position: absolute; 
  top: 165%; 
  left: 50%; 
  transform: translate(-50%);
  background: #FF6700 ;
  padding: 5px;
  border-radius: 3px;
}

</style>
<script setup>
/* eslint-disable */
window.location.href = 'home'
</script>

<template>

  <div class="background"></div>

  <div class="titleWrapper">
      
      <h1 class="titleText">
        <RouterLink to="home"> ADR Prediction</RouterLink>
      </h1>

  </div>

  <RouterLink to="home" style="position: absolute; top: 20%; left: 50%; transform: translate(-50%); color: black;">You will be redirected to the home page shortly</RouterLink>
</template>

<style>

</style>
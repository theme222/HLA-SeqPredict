<script setup>
/* eslint-disable */
import { getAccountInfo, getSequences } from "@/scripts/AccountFunc";
import { ref } from "vue";
import axios from "axios";

const accountName = ref();
const sequenceDropdown = ref([]);
getAccountInfo().then((name) => {
  accountName.value = name;
});

async function addToSequenceList() {
  let serverData = await getSequences();
  for (let sequence of serverData) {
    sequenceDropdown.value.push({ value: sequence[1] });
  }
  console.log(sequenceDropdown.value);
}

addToSequenceList();

</script>

<template>
  <div class="background"></div>

  <div class="titleWrapper">
    <h1 class="titleText">
      <RouterLink to="home">HLA-ADR Prediction</RouterLink>
    </h1>

    <div class="accountName">{{ accountName || "No account" }}</div>
  </div>

  <div class="dashboardHeadingWrapper">
    <h2 class="dashboardHeading">Dashboard</h2>
  </div>

  <div class="dropdownDiv">
    <select name="sequences" id="sequences" class="dropdown">
      <option value="" disabled selected>Select a sequence</option>
      <option v-for="i in sequenceDropdown" :key="i" :value="i.value">
        {{ i.value }}
      </option>
    </select>
    <img
      src="@/assets/down-arrow.svg"
      alt="down-arrow"
      style="
        position: absolute;
        filter: invert(99%) sepia(1%) saturate(594%) hue-rotate(271deg)
          brightness(114%) contrast(100%);
        left: 86%;
        top: 50%;
        width: 25%;
        height: 100%;
        pointer-events: none;
        transform: translate(-50%, -50%);
      "
    />
  </div>

  <button class="igvButton" id="visualizeIGV">View</button>



  <h2 class="tableHeader">Typing results</h2>

  <table class="container">
    <thead>
      <tr>
        <th>Tools</th>
        <th>Gene 1</th>
        <th>ADR 1</th>
        <th>Gene 2</th>
        <th>ADR 2</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td><a href="https://github.com/DiltheyLab/HLA-LA">HLA-LA</a></td>
        <td>A*11:01:01G</td>
        <td>6369</td>
        <td>A*11:01:01G</td>
        <td>01:32:50</td>
      </tr>
      <tr>
        <td><a href="https://github.com/FRED-2/OptiType/">Optitype</a></td>
        <td>9518</td>
        <td>6369</td>
        <td>01:32:50</td>
        <td>01:32:50</td>
      </tr>
      <tr>
        <td><a href="https://github.com/DaehwanKimLab/hisat-genotype">HISAT-genotype</a></td>
        <td>9518</td>
        <td>6369</td>
        <td>01:32:50</td>
        <td>01:32:50</td>
      </tr>
      <tr>
        <td><a href="">HLA Type Bridging</a></td>
        <td>9518</td>
        <td>6369</td>
        <td>01:32:50</td>
        <td>01:32:50</td>
      </tr>
    </tbody>
  </table>

  <div style="position: absolute; top: 100%; height: 7%; width: 20%;"></div>
</template>

<style>

.tableHeader{
  font-family: Arial, Helvetica, sans-serif;
  position: absolute;
  color:black;
  font-size: 2em;
  top: 65%;
  left: 50%;
  transform: translate(-50%);
}

.container {
  font-family: Arial, Helvetica, sans-serif;
  text-align: left;
  overflow: hidden;
  width: 97%;
  display: table;
  height: 29.9%;
  position: absolute;
  background: #36382e;
  color: white;
  left: 50%;
  top: 90%;
  transform: translate(-50%, -50%);
  border-radius: 5px;
}

.container th {
  text-align: center;
  background-color: #6a66a3;
}

/* Background-color of the odd rows */
.container tr:nth-child(odd) {
  background-color: #6d9f71;
}

/* Background-color of the even rows */
.container tr:nth-child(even) {
  background-color: #337357;
}

.container tr td{
  padding-left: 2%;
}

.container tr td a{
  text-decoration: underline ;
  font-weight: 900;
}

.container tr td:hover{
  transition-duration: 0.4s;
  background-color: #36382e;
}


.igvButton {
  position: absolute;
  top: 18.92%;
  left: 300px;
  display: inline-block;
  outline: none;
  cursor: pointer;
  font-weight: 600;
  height: 3.9%;
  border-radius: 3px;
  padding: 12px 24px;
  border: 0;
  color: #3a4149;
  background: #ff6700;
  line-height: 1.15;
  font-size: 1em;
  transform: translate(-50%, -50%);
  width: 90px;
}

.igvButton:hover {
  transition: all 0.1s ease;
  box-shadow: 0 0 0 0 #fff, 0 0 0 3px #1de9b6;
}

.dropdown {
  /* Reset */
  appearance: none;
  border: 0;
  outline: 0;
  font: inherit;
  /* Personalize */
  width: 100%;
  padding: 1rem 4rem 1rem 1rem;
  background: #6a66a3;
  color: white;
  font-family: Arial, Helvetica, sans-serif;
  font-weight: 900;
  border-radius: 0.25em;
  box-shadow: 0 0 1em 0 rgba(0, 0, 0, 0.2);
  cursor: pointer;
}

.dropdownDiv {
  position: absolute;
  top: 17%;
  width: 230px;
}

.igvDivDivlol {
  position: absolute;
  width: 99.4%;
  height: 32%; 
  top:39%;
  left: 50%; 
  transform: translate(-50%,-50%);
  
  background: white;

  box-shadow: 0 0 0 0 #fff, 0 0 0 5px #6a66a3;
  border-radius: 5px;
}

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

.accountName {
  text-wrap: nowrap;
  left: 92%;
  transition: all 0.2s ease;
  position: absolute;
  font-weight: 200;
  color: white;
  font-size: 1.5em;
  top: 17%;
  transform: translate(-50%);
  background: #36382e;
  padding: 0.5%;
  padding-left: 1%;
  padding-right: 1%;
  border-radius: 5px;
}

.accountName:hover {
  transition: all 0.2s ease;
  background: #6a66a3;
  cursor: pointer;
}

.dashboardHeadingWrapper {
  position: absolute;
  top: 6%;
  left: 50%;
}

.dashboardHeading {
  font-family: Arial, Helvetica, sans-serif;
  transform: translate(-50%);
  font-size: 3em;
}

.buttonWrapper {
  position: absolute;
  top: 20%;
  left: 50%;
  transform: translate(-50%);
}

.inputfile {
  width: 0.1px;
  height: 0.1px;
  opacity: 0;
  overflow: hidden;
  position: absolute;
  z-index: -1;
}

.inputfile + label {
  transform: translate(1%);
  padding: 10%;
  padding-left: 20%;
  padding-right: 20%;
  font-size: 1.7em;
  font-weight: 700;
  color: white;
  background-color: #c0c0c0;
  display: inline-block;
  border-radius: 4px;
  border-style: none;
  text-align: center;
}

.inputfile:focus + label,
.inputfile + label {
  cursor: pointer; /* "hand" cursor */
}

.submitButton {
  position: absolute;
  top: 156%;
  left: 50%;
  transform: translate(-50%);
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
}

.submitButton:hover {
  transition: all 0.1s ease;
  box-shadow: 0 0 0 0 #fff, 0 0 0 3px #1de9b6;
}
.filenameDisplay {
  position: absolute;
  left: 50%;
  top: 115%;
  transform: translate(-50%);
  white-space: nowrap;
}

.uploadedFileViewer {
  background-color: rgb(161, 61, 99);
  color: #fff;
  position: absolute;
  padding: 2%;
  border-radius: 4px;
  border-style: none;
  white-space: nowrap;
  top: 15%;
  left: 50%;
  transform: translate(-50%);
  font-size: 2em;
}
</style>
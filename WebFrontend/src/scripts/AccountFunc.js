/* eslint-disable*/
import axios from 'axios'
import Cookies from 'js-cookie'
import {useRouter, useRoute} from 'vue-router'

const router = useRouter()
const route = useRoute()
const backendLink = process.env.VUE_APP_BACKEND_LINK


window.sharedData = {"backendLink":backendLink}
console.log(backendLink)

async function getAccountInfo() {
    if (Cookies.get('session_cookie') == null) {
        console.log("No account found")
        return null
    }
    let response = await axios.post(backendLink+"/api/getAccountInfo", {session_cookie: Cookies.get('session_cookie')}).catch(error => {
        console.error('Error getting info :', error);
        alert("Backend server unavailable")

    });
    return response.data
}

async function getSequences() {
    if (Cookies.get('session_cookie') == null) {
        console.log("No account found")
        return null
    }
    let response = await axios.post(backendLink+"/api/getSequences", {session_cookie: Cookies.get('session_cookie')}).catch(error => {
        console.error('Error getting info :', error);
        alert("Backend server unavailable")
    });
    return response.data['list']
}

export {
    getAccountInfo,
    getSequences,
    backendLink
}

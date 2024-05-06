/* eslint-disable*/
import axios from 'axios'
import Cookies from 'js-cookie'
import {useRouter, useRoute} from 'vue-router'

const router = useRouter()
const route = useRoute()

async function getAccountInfo() {
    if (Cookies.get('session_cookie') == null) {
        console.log("No account found")
        return "No account found"
    }
    let response = await axios.post("http://localhost:7000/api/getAccountInfo", {session_cookie: Cookies.get('session_cookie')}).catch(error => {
        console.error('Error getting info :', error);
        alert("Backend server unavailable")

    });
    return response.data['name']
}

async function getSequences() {
    if (Cookies.get('session_cookie') == null) {
        console.log("No account found")
        return "No account found"
    }
    let response = await axios.post("http://localhost:7000/api/getSequences", {session_cookie: Cookies.get('session_cookie')}).catch(error => {
        console.error('Error getting info :', error);
        alert("Backend server unavailable")
        return 'asdf'
    });
    return response.data['list']
}

export {
    getAccountInfo,
    getSequences
}

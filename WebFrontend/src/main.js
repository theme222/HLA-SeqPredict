import { createApp } from 'vue'
import App from './App.vue'
import router from './router'
import './assets/tailwind.css'

window.addEventListener('error', (event) => {
  window.sharedData.latestError = event.error
  console.error(event.error);
});
const app = createApp(App)
app.use(router).mount('#app')


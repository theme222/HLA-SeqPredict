/* eslint-disable */
import { createApp } from 'vue'
import App from './App.vue'
import * as echarts from 'echarts';
import 'echarts/lib/chart/sunburst';
import router from './router'
import './assets/tailwind.css'

window.addEventListener('error', (event) => {
  window.sharedData.latestError = event.error
  console.error(event.error);
});

const app = createApp(App)
app.config.globalProperties.$echarts = echarts;
app.use(router).mount('#app')


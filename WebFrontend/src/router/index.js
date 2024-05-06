/* eslint-disable */
import { createRouter, createWebHistory } from 'vue-router'
import login from '../views/LoginView.vue'
import home from '../views/HomeView.vue'
import signup from '../views/SignupView.vue'
import upload from '../views/UploadView.vue'
import dashboard from '../views/DashboardView.vue'
import RedirectToHomePlz from '@/scripts/RedirectToHomePlz.vue'

const routes = [
  {path : '/', name: 'redirectToHome', component: RedirectToHomePlz},
  {path : '/home', name: 'home', component: home},
  {path : '/login', name: 'login', component: login},
  {path : '/signup', name: 'signup', component: signup},
  {path : '/upload', name: 'upload', component: upload},
  {path : '/dashboard', name: 'dashboard', component: dashboard}
]

const router = createRouter({
  history: createWebHistory(process.env.BASE_URL),
  routes
})

export default router

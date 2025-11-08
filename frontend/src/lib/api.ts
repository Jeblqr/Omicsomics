import axios from 'axios';

// Use a relative API path so Vite dev server proxy (configured in `vite.config.ts`) can
// forward requests to the backend during development. In other environments a
// deploy-time variable `VITE_API_BASE_PATH` can override the default prefix.
const API_BASE_PATH = (import.meta.env.VITE_API_BASE_PATH as string) || '/api';

// Final axios baseURL should point at the API version prefix. Example: "/api/v1"
const baseURL = `${API_BASE_PATH.replace(/\/$/, '')}/v1`;

console.log(`API requests will be sent to: ${baseURL}`);

const api = axios.create({
  baseURL,
  timeout: 15000,
});

// Add a request interceptor to include auth token and log requests
api.interceptors.request.use(config => {
  // Add auth token from localStorage if available
  const token = localStorage.getItem('auth_token');
  if (token) {
    config.headers.Authorization = `Bearer ${token}`;
  }
  
  console.log(`Requesting URL: ${config.baseURL}${config.url}`);
  return config;
});

// Add a response interceptor to log errors
api.interceptors.response.use(
  response => response,
  error => {
    if (error.config) {
      console.error('API Error:', error.message, '\nURL:', error.config.baseURL + error.config.url);
    } else {
      console.error('API Error:', error);
    }
    if (error.response) {
      console.error('Response:', error.response.status, error.response.data);
    }
    return Promise.reject(error);
  }
);

export default api;

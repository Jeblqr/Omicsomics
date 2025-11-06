import axios from 'axios';

const ROOT_API_URL = (import.meta.env.VITE_API_BASE_URL ?? 'http://localhost:8001').replace(
  /\/$/,
  '',
);

const api = axios.create({
  baseURL: `${ROOT_API_URL}/api/v1`,
  timeout: 15000,
});

export default api;

import { Suspense } from 'react';
import { Navigate, Outlet, Route, Routes } from 'react-router-dom';
import DashboardPage from './dashboard/DashboardPage';
import ProjectsPage from './projects/ProjectsPage';
import RunsPage from './runs/RunsPage';
import DataPage from './data/DataPage';
import SettingsPage from './settings/SettingsPage';
import AuthPage from './auth/AuthPage';
import Layout from '../components/Layout';
import LoadingView from '../components/LoadingView';
import ProtectedRoute from '../components/ProtectedRoute';
import { useAuth } from '../contexts/AuthContext';

const App = () => {
  const { isAuthenticated } = useAuth();

  return (
    <Routes>
      <Route 
        element={isAuthenticated ? <Navigate to="/" replace /> : <AuthPage />} 
        path="/auth" 
      />
      <Route
        element={
          <ProtectedRoute>
            <Layout>
              <Suspense fallback={<LoadingView />}>
                <Outlet />
              </Suspense>
            </Layout>
          </ProtectedRoute>
        }
      >
        <Route element={<DashboardPage />} index />
        <Route element={<ProjectsPage />} path="projects" />
        <Route element={<RunsPage />} path="runs" />
        <Route element={<DataPage />} path="data" />
        <Route element={<SettingsPage />} path="settings" />
      </Route>
    </Routes>
  );
};

export default App;
